from logging import raiseExceptions
import numpy as np
from scipy.io import loadmat
from scipy.fft import fft, ifft

def covariance_from_array(array, periodic=False):
    halflen = int(np.floor(np.shape(array)[0]/2))
    if periodic: halflen = int(0)
    array = array - array[halflen::].mean()
    rms = np.std(array[halflen::], ddof = 1)
    err = []
    # rms = np.sqrt(np.sum(data[halflen::]**2)/(halflen-1))
    for step in range(len(array)- halflen):
        cov = array * np.roll(array, step) /rms/rms
        err.append(np.average(cov[halflen::]))
    return err

def covariance_2d(data_2d, periodic=False):
    err = 0
    for ind, array in enumerate(data_2d):
        err += np.array( covariance_from_array(array, periodic))

    err = err/(ind+1)
    return err
    
    
# def params_finite_difference(M:int, x0:float, x:np.array):
#     """_summary_

#     Args:
#         M (int): order of derivative.
#         x0 (float): _description_
#         x (np.array): _description_
#     """
#     N = x.shape[0]-1
    
#     params = np.empty(shape=(M+1, N+1, N+1)) # indices (m, n, nu)
#     params[:] = np.nan
    
#     params[0, 0, 0] = 1
#     c1 = 1
#     for n in range(1, N+1):
#         c2 = 1
#         for nu in range(n):
#             c3 = x[n] - x[nu]
#             c2 = c2 * c3
#             if n<M: params[n, n-1, nu] = 0
#             for m in range(min(n, M)+1):
#                 params[m, n, nu] = ((x[n] - x0) * params[m, n-1, nu] - m * params[m-1, n-1, nu]) / c3
        
#         for m in range(min(n, M)+1):
#             params[m, n, n] = c1/c2 * (m*params[m-1, n-1, n-1] - (x[n-1] - x0) * params[m, n-1, n-1])
        
#         c1 = c2
    
#     return params
    
    

class specturm():
    """_summary_
    """
    def __init__(self, data: np.ndarray):
        self.data = data
        self.ndims = data.ndim()
        self.shape = data.shape()
        return 
    
        
    
if __name__ == "__main__":
    x = np.arange(-4, 4)
    x0 = 0
    M = 4
    params = params_finite_difference(M, x0, x)
    print(params)