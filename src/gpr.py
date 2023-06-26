from logging import raiseExceptions
import numpy as np
from scipy.io import loadmat
from scipy.fft import fft, ifft

ODI_path = '/Users/user/Documents/Projects/GaussianProcessRegression/JHTDB_GPR/outputs/Realizations_tecplot/ODI/data%.4i_ODI.mat'

class radial():
    """
    Gaussian Process Regression with Radial Basis Kernel
    """
    def __init__(self, ):
        self.Kxx       = lambda x1, y1, x2, y2, l: (1/l**2 - (x1[:, np.newaxis] - x2[np.newaxis, :])**2/(l**4)) * np.exp(-0.5*((x1[:, np.newaxis] - x2[np.newaxis, :])**2+(y1[:, np.newaxis] - y2[np.newaxis, :])**2)/(l**2))
        self.Kyy       = lambda x1, y1, x2, y2, l: (1/l**2 - (y1[:, np.newaxis] - y2[np.newaxis, :])**2/(l**4)) * np.exp(-0.5*((x1[:, np.newaxis] - x2[np.newaxis, :])**2+(y1[:, np.newaxis] - y2[np.newaxis, :])**2)/(l**2))
        self.Kxy       = lambda x1, y1, x2, y2, l: -((x1[:, np.newaxis] - x2[np.newaxis, :]) * (y1[:, np.newaxis] - y2[np.newaxis, :]))/l**4 * np.exp(-0.5*((x1[:, np.newaxis] - x2[np.newaxis, :])**2+(y1[:, np.newaxis] - y2[np.newaxis, :])**2)/(l**2))
        self.Kx2       = lambda x1, y1, x2, y2, l: (x1[:, np.newaxis] - x2[np.newaxis, :]) * np.exp(-0.5*((x1[:, np.newaxis] - x2[np.newaxis, :])**2 + (y1[:, np.newaxis] - y2[np.newaxis, :])**2)/(l**2))
        self.Ky2       = lambda x1, y1, x2, y2, l: (y1[:, np.newaxis] - y2[np.newaxis, :])/(l**2) * np.exp(-0.5*((x1[:, np.newaxis] - x2[np.newaxis, :])**2 + (y1[:, np.newaxis] - y2[np.newaxis, :])**2)/(l**2))
        self.K         = lambda x1, y1, x2, y2, l: np.exp(-0.5*((x1[:, np.newaxis] - x2[np.newaxis, :])**2 + (y1[:, np.newaxis] - y2[np.newaxis, :])**2)/(l**2))


        self.Krr       = lambda x1, y1, x2, y2, l: (1/l**2 - ((x1[:, np.newaxis] - x2[np.newaxis, :])**2+(y1[:, np.newaxis] - y2[np.newaxis, :])**2)/(l**4)) * np.exp(-0.5*((x1[:, np.newaxis] - x2[np.newaxis, :])**2+(y1[:, np.newaxis] - y2[np.newaxis, :])**2)/(l**2))

    def radial_kernel(self, r1, r2, l):
        """
        radial basis kernel k = exp(-|r|^2/(2l^2)

        input: 
            r1: ndarray(N, d) N position vectors of d dimensions
            r2: ndarray(M, d) M position vectors of d deminsions
            l: length scale

        return:
            k: ndarray(N*d, M*d) covariance matrix of r1 and r2
        """
        


        k = np.exp(0)
        return k
    
    def regression(self, Xstar, Ystar, X, Y, dPdx, dPdy, l, noise):
        observation = np.concatenate((dPdx, dPdy, np.array([0])))
        
        s11 = self.Kxx(X, Y, X, Y, l) + noise**2*np.eye(np.shape(X)[0])
        s12 = self.Kxy(X, Y, X, Y, l)
        s13 = -self.Kx2(X, Y, np.array([0]), np.array([0]), l)

        s21 = self.Kxy(X, Y, X, Y, l)
        s22 = self.Kyy(X, Y, X, Y, l) + noise**2*np.eye(np.shape(X)[0])
        s23 = -self.Ky2(X, Y, np.array([0]), np.array([0]), l)


        s31 = self.Kx2(np.array([0]), np.array([0]), X, Y, l)
        s32 = self.Ky2(np.array([0]), np.array([0]), X, Y, l)
        s33 = 1

        sigma11 = np.block([[s11, s12, s13 ], [s21, s22, s23], [s31, s32, s33]])
        sigma21 = np.block([self.Kx2(Xstar, Ystar, X, Y, l), self.Ky2(Xstar, Ystar, X, Y, l), self.K(Xstar, Ystar, np.array([0]), np.array([0]), l)])
        sigma12 = np.block([[-self.Kx2(X, Y, Xstar, Ystar, l)], [-self.Ky2(X, Y, Xstar, Ystar, l)], [self.K(np.array([0]), np.array([0]), Xstar, Ystar,  l)]])

        k_p = self.K(Xstar, Ystar, Xstar, Ystar, l)

        mu =  sigma21 @ np.linalg.solve(sigma11, observation)
        std = k_p - sigma21 @ np.linalg.solve(sigma11, sigma12)
        return mu, std

class GPR():
    """
    class of Gaussian Process Regression
    """
    def __init__(self, dims=150):
        """
        dims: dimension of GPRP
        """
        # self.ODIP_path = '/Users/user/Documents/Projects/GaussianProcessRegression/JHTDB_GPR/outputs/Realizations_tecplot/ODI/data0001_ODI.mat'
        # self.GPRP_path = '/Users/user/Documents/Projects/GaussianProcessRegression/JHTDB_GPR/outputs/results_tecplot/Variance40/Realizations_tecplot/Optimized_l6.2500e-02_data0001.mat'
        self.REALP_path = '/Users/user/Documents/Projects/GaussianProcessRegression/JHTDB_GPR/inputs/realP.mat'
        self.XY_path = '/Users/user/Documents/Projects/GaussianProcessRegression/JHTDB_GPR/inputs/XY.mat'

        self.dims = int(dims)

        self.realP = self.readdata( self.REALP_path)
        self.stdP = np.std(self.realP, ddof = 1 )

        # read coordinate system
        coord = loadmat(self.XY_path)
        self.X = coord['X_samples'][:self.dims, :self.dims]
        self.Y = coord['Y_samples'][:self.dims, :self.dims]
        # move to zero
        self.X -= np.min(self.X)
        self.Y -= np.min(self.Y)

        self.dx = self.X[1, 0] - self.X[0, 0]

        # Kronmogorov Length Scale
        self.kdelta = 0.00280
        self.transpose = True

    def readdata(self, path, transpose=False):
        P = readmat(path)
        ddims = int(np.sqrt(len(P.flatten())))
        P = np.reshape(P, (ddims, ddims))[:self.dims, :self.dims]
        if transpose: P = P.T
        self.P = P - np.mean(P)
        return self.P



    def read_GPR(self, path, realization=1):
        self.GPR_path = path
        self.GPRP = self.readdata(path % realization, self.transpose)
        return self.GPRP

    def read_ODI(self, realization=1, path=ODI_path, ):
        self.ODI_path = path
        self.ODIP = self.readdata(path % realization, self.transpose)
        return self.ODIP
    
    def compare(self, realization=1):
        self.read_GPR(self.GPR_path, realization)
        self.read_ODI(path = self.ODI_path, realization=realization)
        print('load realization %.4i to compare' % realization)
        return self

    def error_normalized(self, P):
        n_err = np.std(P-self.realP, ddof=1)/self.stdP
        return n_err