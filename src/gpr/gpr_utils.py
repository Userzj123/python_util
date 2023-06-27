from scipy.optimize import curve_fit




def fit_length_scale_pressure(x, data, kernel):
    popt, pcov = curve_fit(kernel, x, data, bounds=(0, 1))
    l:float = popt
    
    return l 