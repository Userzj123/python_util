import numpy as np 
import matplotlib.pyplot as plt
from .generate_field import *
from matplotlib import cm, colors
import os
import imageio


def contour_channel(domain, dims, data, colormap = cm.Oranges):
    """
    plot the contour of x, y, z cross section in channel flow at the maximum index.

    Args:
        domain (_type_): _description_
        dims (_type_): _description_
        data (_type_): _description_

    Returns:
        _type_: _description_
    """
    from matplotlib.gridspec import GridSpec
    import matplotlib.ticker as tick
    
    x_coords, y_coords, z_coords = xyz(domain, dims)
    x_ind, y_ind, z_ind = np.unravel_index(data.argmax(), dims)
    
    vmin = round(data.min(), 1)
    vmax = round(data.max(), 1)
    levels = np.linspace(vmin, vmax, 101)

    fig = plt.figure(layout="constrained", figsize=(8,4))

    gs = GridSpec(2, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    # identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
    ax2 = fig.add_subplot(gs[1, :])
    ax3 = fig.add_subplot(gs[0, 1])
    axes = [ax1, ax2, ax3]


    images = []
    images.append(ax1.contourf(x_coords, y_coords, data[:, :, z_ind].T, levels, cmap=colormap, extend='both'))
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_aspect('equal', 'box')
    ax1.set_ylim(bottom=0)
    # cbar = fig.colorbar(cs)


    images.append(ax2.contourf(x_coords, z_coords, data[:, y_ind, :].T, levels, cmap=colormap, extend='both'))
    ax2.set_xlabel('x')
    ax2.set_ylabel('z')
    ax2.set_aspect('equal', 'box')
    ax2.set_ylim(bottom=0)
    # cbar = fig.colorbar(cs1)

    images.append(ax3.contourf(y_coords, z_coords, data[x_ind, :, :].T, levels, cmap=colormap, extend='both'))
    ax3.set_xlabel('y')
    ax3.set_ylabel('z')
    ax3.set_aspect('equal', 'box')
    ax3.set_ylim(bottom=0)
    # cbar = fig.colorbar(cs2)


    cax = plt.axes([1.05, 0.11, 0.025, 0.78])
    cbar = fig.colorbar(images[1], cax=cax)
    # cbar.set_ticks(levels[::10])
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))

    fig.show()

    return fig, axes