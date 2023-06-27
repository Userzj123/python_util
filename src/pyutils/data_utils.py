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
    
    
class specturm():
    """_summary_
    """
    def __init__(self, data: np.ndarray):
        self.data = data
        self.ndims = data.ndim()
        self.shape = data.shape()
        return 
    
        
    
if __name__ == "__main__":
    xx = np.linspace(0, 2*np.pi, 128)
    f = lambda x: np.sin(x)
    covariance_from_array(f(xx))