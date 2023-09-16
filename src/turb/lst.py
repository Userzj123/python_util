
import numpy as np
from scipy import linalg


class lst():
    def __init__(self, N:int, bf:int) -> None:
        """Initialize the LST Analysis

        Args:
            N (int): number of Chebyshev polynomials
            bf (int): set = 1 for Couette, 2 for Poiseuille, 3 for quiescent
            
        Returns:

        """
        # Input parameters 
        self.N = N          # number of Chebyshev polynomials

        R  = 180         # Reynolds number
        kx = 1           # streamwise wavenumber
        kz = 0           # spanwise wavenumber
        Ri = 0.0        # Richardson number
        Pr = 0.71        # Prantl number
        Ra = 8*R*R/Pr*Ri # Rayleigh number
        print('Rayleigh number is '+ str(Ra))

        bf = 2   # 
        
        
        # Set up grid and differentiation matrices
        self.y_phys                = np.cos(np.linspace(0, np.pi, N))[:, np.newaxis]   # Generate Chebyshev grid for base flow solver
        self.D0, self.D1, self.D2, self.D3, self.D4 = self.dmat(N)   # Chebyshev polynomials and derivatives at the Gauss points
        
        # Find the base flow
        [U,Up,Upp,T,Tp] = self.bounded_base(self.y_phys,N,bf)


        # Find eigenvalues of stability operators
        self.A, self.B = self.Operator(kx,kz,R,Pr,Ri,U,Up,Upp,Tp,self.D0, self.D1, self.D2, self.D4)

        self.omega, self.q = self.solve_eig(self.A, self.B)

        # [WorkingOn]
        # # Find fields in physical space
        # [u,v,w,eta] = vel_field(kx,kz,D0p,D1p,q,N)
        # Temperature = D0p*q(2*N+1:3*N,:)
        return
    
    def solve_eig(self, A, B):
        # find eigenvalues
        # omega, q = linalg.eig(A, B)
        omega, q = np.linalg.eig(np.linalg.inv(B)@A)

        omega = 1j*omega   # eigenvalues omega in vector form

        # remove bad eigenvalues
        sp = np.logical_and(abs(omega)>1e-10, abs(omega)<50)

        omega = omega[sp]
        q = q[:, sp]
        return omega, q
    
    
    def dmat(self, N):
        num = N-1
        D0 = np.cos(np.arange(N)[np.newaxis, :]  * np.pi * np.arange(N)[:, np.newaxis] / num )
        
        # create higher derivative matrices
        D1 = np.concatenate((np.zeros(shape=(N, 1)), D0[:, 0][:, np.newaxis], 4*D0[:, 1][:, np.newaxis]), axis=1)
        D2 = np.concatenate((np.zeros(shape=(N, 2)),                          4*D1[:, 1][:, np.newaxis]), axis=1)
        D3 = np.zeros(shape=(N, 3))
        D4 = np.zeros(shape=(N, 3))



        for j in range(3, N):
            D1= np.concatenate((D1, 2*j*D0[:, j-1][:, np.newaxis]+j*D1[:, j-2][:, np.newaxis]/(j-2)), axis=1)
            D2= np.concatenate((D2, 2*j*D1[:, j-1][:, np.newaxis]+j*D2[:, j-2][:, np.newaxis]/(j-2)), axis=1)
            D3= np.concatenate((D3, 2*j*D2[:, j-1][:, np.newaxis]+j*D3[:, j-2][:, np.newaxis]/(j-2)), axis=1)
            D4= np.concatenate((D4, 2*j*D3[:, j-1][:, np.newaxis]+j*D4[:, j-2][:, np.newaxis]/(j-2)), axis=1)
            
        return D0, D1, D2, D3, D4
    
    def bounded_base(self, y_phys,N,bf):
        if bf == 1: # Couette flow
    
            U   = y_phys
            Up  = np.ones(shape=(N, 1))
            Upp = np.zeros(shape=(N, 1))
            T   = y_phys
            Tp  = np.zeros(shape=(N, 1))
            
        elif bf == 2 : # Poiseuille flow
            
            U   = (1 - y_phys**2)
            Up  = -2*y_phys
            Upp = -2*np.ones(shape=(N, 1))
            T   = np.zeros(shape=(N, 1))
            Tp  = np.zeros(shape=(N, 1))
            
        elif bf == 3: # quiescent flow
            
            U   = np.zeros(shape=(N, 1))
            Up  = np.zeros(shape=(N, 1))
            Upp = np.zeros(shape=(N, 1))
            T   = y_phys
            Tp  = np.ones(shape=(N, 1))
        else:
            print('Need to select bf = 1 or 2')
            
            
        return U,Up,Upp,T,Tp
    
    def Operator(self, kx,kz,R,Pr,Ri,U,Up,Upp,Tp,D0,D1,D2,D4):
        # Useful variables --------------------------------------------------------
        k2 = kx**2 + kz**2       # wavenumber^2
        N  = U.shape[0]
        M  =  np.ones(shape=(1, N)) # matrix for mean flow variables 
        er = -200*1j           # for spurious eigenvalues from BCs

        LSQ = -1j*kx*(U@M)*D0 + (1/R)*(D2-k2*D0)
        LOS = -1j*kx*(U@M)*(D2-k2*D0) + 1j*kx*(Upp@M)*D0 + (1/R)*(D4-(2*k2*D2)+((k2**2)*D0))
        A32 = -1j*kx*(U@M)*D0 + (1/R/Pr)*(D2-k2*D0)
        
        A = np.block([[LOS, np.zeros((N, N)), -Ri*k2*D0],[-1j*kz*(Up@M*D0) , LSQ, np.zeros((N, N))],[-(Tp@M*D0), np.zeros((N, N)), A32]])

        B = np.block([[D2-k2*D0, np.zeros((N, N)), np.zeros((N, N))], [np.zeros((N, N)), D0, np.zeros((N, N))], [np.zeros((N, N)),  np.zeros((N, N)), D0]])

        B[0, 0:N] = D0[0, :]
        B[1, 0:N, ] = D1[0, :] # v, Dv at y = top
        B[N-2, 0:N] = D1[N-1, :]
        B[N-1, 0:N] = D0[N-1, :] # v, Dv at y = bot

        # apply boundary condtions to top and bottom 2 rows (i.e. v = Dv = 0)
        A[0, 0:N] = er*D0[0, :]       # vanishing in the free stream
        A[1, 0:N] = er*D1[0, :]       # vanishing in the free stream
        A[N-2, 0:N] = er*D1[N-1, :]   # gradient vanishing at wall (no slip)
        A[N-1, 0:N] = er*D0[N-1, :]   # no penetration at wall

        # clear the rows to apply Squire boundary conditions ----------------------
        A[N, :] = 0
        A[2*N-1, :] = 0

        # apply Squire boundary conditions (eta = 0 at y=top,bot)
        A[N, N:2*N, ] = er*D0[0, :]
        A[2*N-1, N:2*N, ] = er*D0[N-1, :]

        # apply temperature boundary conditions
        A[2*N, :] = 0
        A[3*N-1, :] = 0

        A[2*N, 2*N:3*N] = er*D0[0, :]
        A[3*N-1, 2*N:3*N] = er*D0[N-1, :]
        return A, B
    
    
if __name__ == "__main__":
    channel = lst(384, 2)
    
    print('test')