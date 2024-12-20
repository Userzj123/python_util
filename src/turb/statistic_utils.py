import numpy as np
import matplotlib.pyplot as plt

def corr_multidimensional(arrays, periodic=False, num_step=-1):
    """ 
        From arrays with dims = [..., ind_sampling] to obtain correlation with dims [..., ind_d ]
    """
    origin_shape = arrays.shape

    # Determine sampling length
    halflen = int(np.floor(arrays.shape[-1]/2))
    if periodic: halflen = int(0)

    # Resize arrays
    arrays = arrays.reshape([-1, origin_shape[-1]])

    mean = arrays[:, halflen:].mean(axis=-1)[:, np.newaxis]
    arrays = arrays - mean
    rms = np.std(arrays[:, halflen:], ddof = 0, axis=-1)[:, np.newaxis]

    if isinstance(num_step, int) and num_step == -1:
        num_step = arrays.shape[-1] - halflen
        
        err = np.zeros((*arrays.shape[:-1], num_step))
        # print(err.shape)
        for step in range(num_step):
            err[:, step] = np.average((arrays * np.roll(arrays, step, axis=-1) /rms/rms)[:, halflen:], axis=-1)

        corr = err.reshape((*origin_shape[:-1], num_step))
    elif isinstance(num_step, int) and num_step > -1:
        err = np.zeros((*arrays.shape[:-1], num_step))
        # print(err.shape)
        for step in range(num_step):
            err[:, step] = np.average((arrays * np.roll(arrays, step, axis=-1) /rms/rms)[:, halflen:], axis=-1)
        corr = err.reshape((*origin_shape[:-1], num_step))
    elif isinstance(num_step, np.ndarray) or isinstance(num_step, list):
        err = np.zeros((*arrays.shape[:-1], len(num_step)))
        # print(err.shape)
        i = 0
        for step in num_step:
            err[:, i] = np.average((arrays * np.roll(arrays, step, axis=-1) /rms/rms)[:, halflen:], axis=-1)
            i += 1
        corr = err.reshape((*origin_shape[:-1], len(num_step)))
    else:
        raise Exception("None of the num_step type is matched.")
    return corr


def plot_symmetric(ax, x, y, **kwargs):
    xx = np.concatenate((-x[::-1], x))
    yy = np.concatenate((y[::-1], y))
    return ax.plot(xx, yy, **kwargs)


def tau_integral(corr, tau):
    """_summary_

    Args:
        corr (ndarray([*, Nt])): Temporal autocorrelation function.
        tau (ndarray([Nt])): Time lag of measurement.

    Returns:
        tau_int: Integral timescale.
    """
    dt = tau[1]- tau[0]
    inds = np.where(corr<0)[0]
    if inds.size == 0 :
        tau_int = np.sum(corr[:]*dt) + np.sum(corr[1:]*dt)
    else:
        tau_int = np.sum(corr[:inds[0]]*dt) + np.sum(corr[1:inds[0]]*dt) 
    return tau_int

def tau_integrals(corrs, t):
    original_shape = corrs.shape

    corrs = corrs.reshape((-1, original_shape[-1]))

    taus = np.zeros((corrs.shape[0]))

    for ind, tau in enumerate(taus):
        taus[ind] = tau_integral(corrs[ind, :], t)

    return taus.reshape(original_shape[:-1])

    
def time_averaged_measure(m, window):
    """
    m: [Nsample, Nmeasurement]
    window: int 
    """

    import math
    
    max_n = math.floor(m.shape[-1]/window)
    # print(max_n, window)
    var_sample = m[:, :max_n*window].reshape((m.shape[0], max_n, -1)).mean(axis=-1)
    return var_sample
