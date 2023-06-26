import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors


def contour_single(x, y, data, levs=61, eshrink=5/6, **fig_kw):
    """
    Generate single contour given x, y and data.
    """
    emax = np.ceil(np.max([data.max(),]))
    emin = np.floor(np.min([data.min(),]))
    levels = np.linspace(emin, emax, levs)

    fig, ax = plt.subplots(figsize=(5,5), dpi = 150,  sharex=True, sharey=True, **fig_kw)

    cs1 = ax.contourf(x, y, data, levels, cmap=cm.twilight_shifted, extend='both')
    ax.set_xlabel(r'$x$', size=15)
    ax.set_ylabel(r'$y$', size=15)

    ax.set_aspect('equal')

    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.825, 0.16, 0.025, 0.68])
    cbar = fig.colorbar(cs1, cax=cax, )
    cbar.set_ticks(levels[::15])
    return fig, ax