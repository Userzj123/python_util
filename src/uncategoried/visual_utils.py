import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.io import loadmat
import math
from matplotlib import cm, colors
from data_utils import *
from scipy.fft import fft, ifft
from scipy.optimize import curve_fit

# font_family = 'Times New Roman'
def tecplot_reader(file, nb_var):
    """Tecplot reader."""
    arrays = []
    with open(file, 'r') as a:
        for idx, line in enumerate(a.readlines()):
            if idx < 3:
                continue
            else:
                try: arrays.append([float(s) for s in line.split()])
                except: print(line)

    arrays = np.concatenate(arrays)
    output = np.split(arrays, nb_var)

    return output


def cumulative_error(x, y, Aerr):
    fig, ax = plt.subplots()
    plt.rc('text', usetex=True)
    # plt.rc('font', family=font_family)
    X, Y = np.meshgrid(x, y)

    lev_exp = np.linspace(np.log10(Aerr.min())-0.1,
                    np.log10(Aerr.max())+0.1, 30)
    levs = np.power(10, lev_exp)
    cs = ax.contourf(X, Y, Aerr, levs, norm=colors.LogNorm(), cmap=cm.bone)
    cbar = fig.colorbar(cs)
    cbar.set_ticks([1, ])
    return fig, ax, cbar


def comparison(gpr, eshrink=1, levP = 11, fontsize = 16, levstep=5):
    real_P = gpr.realP
    GPRP = gpr.GPRP
    ODIP = gpr.ODIP

    stdP = np.std(real_P, ddof = 1 )

    vmin = math.floor(np.min([real_P.min(), GPRP.min(), ODIP.min()]))
    vmax = math.ceil(np.max([real_P.max(), GPRP.max(), ODIP.max()]))
    levels = np.linspace(-1.2*eshrink, 1.2*eshrink, levP)
    # levels = np.linspace(vmin, vmax, levP)

    diffGPR = (GPRP - real_P)/stdP
    diffODI = (ODIP - real_P)/stdP
    diff = (GPRP - ODIP)/stdP

    emax = math.ceil(np.max([np.abs(diffGPR).max(), np.abs(diffODI).max(), np.abs(diff).max()]))
    # emin = math.floor(np.min([diffGPR.min(), diffODI.min(), diff.min()]))
    # emax = math.ceil(np.max([diffGPR.max(), diffODI.max(), diff.max()]))
    elevs = np.linspace(-emax, emax, 21)

    # Plot
    plt.close('all')
    fig, axes = plt.subplots(2, 3, sharex=True, sharey=True)
    axes = axes.flat

    plt.rc('text', usetex=True)
    # plt.rc('font', family=font_family)

    cs1 = axes[0].contourf(gpr.X, gpr.Y, real_P, levels, cmap=cm.twilight_shifted, extend='both')
    cs2 = axes[1].contourf(gpr.X, gpr.Y, GPRP, levels, cmap=cm.twilight_shifted, extend='both')
    cs3 = axes[2].contourf(gpr.X, gpr.Y, ODIP, levels, cmap=cm.twilight_shifted, extend='both')

    cs4 = axes[3].contourf(gpr.X, gpr.Y, diff, elevs, cmap=cm.coolwarm,)
    cs5 = axes[4].contourf(gpr.X, gpr.Y, diffGPR, elevs, cmap=cm.coolwarm,)
    cs6 = axes[5].contourf(gpr.X, gpr.Y, diffODI, elevs, cmap=cm.coolwarm,)

    cbar = fig.colorbar(cs3, ax=[axes[0], axes[1], axes[2]])
    cbar.set_ticks(levels[::levstep])
    cbar.set_label(r'Pressure', size=fontsize)

    ecbar = fig.colorbar(cs6, ax= [axes[3], axes[4], axes[5]],)
    ecbar.set_ticks(elevs[::5])
    ecbar.set_label(r'Error $\epsilon$', size=fontsize)

    axes[0].set_title('a)')
    axes[1].set_title('b)')
    axes[2].set_title('c)')
    axes[3].set_title('d)')
    axes[4].set_title('e)')
    axes[5].set_title('f)')

    # axes[0].set_title('a) True Pressure')
    # axes[1].set_title('b) $P_{GPR}$')
    # axes[2].set_title('c) $P_{ODI}$')
    # axes[3].set_title('d) $P_{GPR}-P_{ODI}$')
    # axes[4].set_title('e) GPR error')
    # axes[5].set_title('f) ODI error')

    axes[0].set_ylabel(r'$y$', fontsize = fontsize)
    axes[3].set_ylabel(r'$y$', fontsize = fontsize)
    axes[0].set_xlabel(r'$x$', fontsize = fontsize)
    axes[1].set_xlabel(r'$x$', fontsize = fontsize)
    axes[2].set_xlabel(r'$x$', fontsize = fontsize)
    axes[3].set_xlabel(r'$x$', fontsize = fontsize)
    axes[4].set_xlabel(r'$x$', fontsize = fontsize)
    axes[5].set_xlabel(r'$x$', fontsize = fontsize)

    # fig.text(0.5, 0.05, 'X', ha='center', size=22)
    # fig.text(0.1, 0.5, 'Y', va='center', rotation='vertical', size=22)
    return fig, axes


def correlation(P, gpr, ax=plt.subplots(), fitting=False, func=radial_basis):
    dx = gpr.dx
    err = []
    for ind, data in enumerate(P):
        halflen = int(np.shape(data)[0]/2)
        data = data - data.mean()
        rms = np.std(data, ddof = 1)
        # rms = np.sqrt(np.sum(data[halflen::]**2)/(halflen-1))
        for step in range(halflen):
            corr = data * np.roll(data, step) /rms/rms
            # err.append(np.average(corr[halflen::]))
            err.append(np.average(corr))

    plt.rc('text', usetex=True)
    # plt.rc('font', family=font_family)

    cov = np.average(np.reshape(err, (ind+1, step+1)), axis=0)
    x = range(len(cov)) * dx


    popt = 0
    if fitting:
        ax.plot(x, cov, label=r'correlation function of $\tilde{p}$', color='black')
        popt, pcov = curve_fit(func, x, cov, bounds=(0, 1))
        print(popt)
        sigma = popt

        ax.plot(x, func(x, sigma), label='fitted correlation function', linestyle='dashdot', color='indianred')

    return cov, err, popt




def plot_energy_spectrum(P, dx, kdelta, ax, label='Real Pressure', color='r'):
    w, s = energy_spectrum(P, dx)

    ax.loglog(w*kdelta, s, label=label, color=color)
    return w, s

if __name__ == "__main__":
    gpr = GPR()
    fig, ax = plt.subplots()
    cov, err, popt = correlation(gpr.realP, gpr, ax, fitting=True)
    



        