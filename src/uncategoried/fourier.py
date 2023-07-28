import numpy as np


def covariance_from_array(data):
        halflen = math.floor(np.shape(data)[0]/2)
        data = data - data.mean()
        rms = np.std(data, ddof = 1)
        # rms = np.sqrt(np.sum(data[halflen::]**2)/(halflen-1))
        for step in range(halflen):
            corr = data * np.roll(data, step) /rms/rms
            # err.append(np.average(corr[halflen::]))
            err = np.average(corr)
        return err

def covariance_2d(data, dx):
    err = []
    for ind, data in enumerate(data):
        err += covariance_from_array(data)

    err = err/(ind+1)
    x = range(len(cov)) * dx
    return err

def power_spectrum_from_array(data, dx):
    N = len(data)
    X = N * dx
    Kx = 2*np.pi/(2*dx)
    DKx = Kx / (N/2)
    half = int(N/2 + 1)
    AngularX = np.linspace(0, (N/2)*DKx, half)

    sigmad2 = np.std(data, ddof=1)**2
    P2 = fft(data)*dx

    P1 = P2[:half]
    S2 = (2/X) *(np.abs(P1)**2)
    return AngularX, S2

def power_spectrum_2d(P, dx):
    w = 0
    s = 0
    for row, data in enumerate(P):
        tmp_w, tmp_s = power_spectrum_from_array(data, dx)
        w += tmp_w
        s += tmp_s
    row += 1
    
    w = w/row
    s = s/row
    return w, s


if __name__ == "__main__":
    test = np.sin(np.linspace(0, 2*np.pi, 100))
    power_spectrum_from_array(test)