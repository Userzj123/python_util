import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors


def generate_gif(timestep, ldata, fig_dir, var, k, marker, **fig_args):
    import os
    import imageio
    gif_dir = './gif'
    if not os.path.exists(gif_dir):
        os.makedirs(gif_dir)
    filenames = []
    for i in timestep:
        ldata.read_data(i)
        fig, ax = contour_channel(ldata.coords, ldata.data[var][k], **fig_args)
        
        # create file name and append it to a list
        filename = gif_dir + f'/%.5i.png' % i
        filenames.append(filename)
        
        # save frame
        fig.savefig(filename, bbox_inches='tight')
        plt.close()
    # build gif
    gif_fname = fig_dir + '/%s_%.2i_%s.gif' %(var, k, marker)
    with imageio.get_writer(gif_fname, mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)
            
    # Remove files
    for filename in set(filenames):
        os.remove(filename)
    return 

def general_gif(tt, plt_func, pltfunc_args, **kwargs):
    import os
    import imageio
    
    defaultKwargs = {
        'gif_fname' : './result.gif' 
    }

    kwargs = { **defaultKwargs, **kwargs }
    
    
    gif_dir = './tmp_gif'
    if not os.path.exists(gif_dir):
        os.makedirs(gif_dir)
    filenames = []
    for it in tt:
        pltfunc_args['lesgo_data'].read_data(it)
        pltfunc_args['lesgo_data'].current_tstep = it
        fig, ax = plt_func(**pltfunc_args)
        
        # create file name and append it to a list
        filename = gif_dir + f'/%.5i.png' % it
        filenames.append(filename)
        
        # save frame
        fig.savefig(filename, bbox_inches='tight')
        plt.close()
    # build gif
    with imageio.get_writer(kwargs['gif_fname'], mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)
            
    # Remove files
    for filename in set(filenames):
        os.remove(filename)
    return 


def contour_single(coords, data, levs=61, eshrink=5/6, **fig_kw):
    """
    Generate single contour given x, y and data.
    """
    x, y = coords
    
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

def contour_channel(coords, data, **kwargs):
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
    
        
    defaultKwargs = {
        'alpha' : 0.1,
        'cmap' : cm.pink_r,
        'tick_fmt' : '%.1f',
        'figsize' : (8, 8),
        'levels': 101,
        'clip'  : False,
        'norm'  : colors.Normalize(data.min(), data.max()),
        'cbar_label' : 'theta',
        'cross_indices' : None
    }
    
    
    kwargs = { **defaultKwargs, **kwargs }
    
    
    x_coords, y_coords, z_coords = coords
    domain = tuple([d.max() for d in coords])
    dims = tuple([d.shape[0] for d in coords])
    if kwargs["cross_indices"] is None:
        x_ind, y_ind, z_ind = np.unravel_index(data.argmax(), dims)
    else:
        x_ind, y_ind, z_ind = kwargs["cross_indices"]
    print("maximum index = [%i, %i, %i]" %(x_ind, y_ind, z_ind))

    fig = plt.figure(layout="constrained", figsize=(8,4))

    gs = GridSpec(2, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    # identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
    ax2 = fig.add_subplot(gs[1, :])
    ax3 = fig.add_subplot(gs[0, 1])
    axes = [ax1, ax2, ax3]


    images = []
    ax1.vlines(x_coords[x_ind], y_coords.min(), y_coords.max(), color='navy', linestyles='dashdot', alpha=kwargs['alpha'])
    ax1.hlines(y_coords[y_ind], x_coords.min(), x_coords.max(), color='navy', linestyles='dashdot', alpha=kwargs['alpha'])
    images.append(ax1.contourf(x_coords, y_coords, data[:, :, z_ind].T, levels=kwargs['levels'],
                    norm=kwargs['norm'], cmap=kwargs['cmap'], extend='both'))
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_aspect('equal', 'box')
    ax1.set_ylim(bottom=0)
    # cbar = fig.colorbar(cs)

    ax2.vlines(x_coords[x_ind], z_coords.min(), z_coords.max(), color='navy', linestyles='dashdot', alpha=kwargs['alpha'])
    ax2.hlines(z_coords[z_ind], x_coords.min(), x_coords.max(), color='navy', linestyles='dashdot', alpha=kwargs['alpha'])
    images.append(ax2.contourf(x_coords, z_coords, data[:, y_ind, :].T, levels=kwargs['levels'],
                    norm=kwargs['norm'], cmap=kwargs['cmap'], extend='both'))
    ax2.set_xlabel('x')
    ax2.set_ylabel('z')
    ax2.set_aspect('equal', 'box')
    ax2.set_ylim(bottom=0)
    # cbar = fig.colorbar(cs1)

    ax3.vlines(y_coords[y_ind], z_coords.min(), z_coords.max(), color='navy', linestyles='dashdot', alpha=kwargs['alpha'])
    ax3.hlines(z_coords[z_ind], y_coords.min(), y_coords.max(), color='navy', linestyles='dashdot', alpha=kwargs['alpha'])
    images.append(ax3.contourf(y_coords, z_coords, data[x_ind, :, :].T, levels=kwargs['levels'],
                    norm=kwargs['norm'], cmap=kwargs['cmap'], extend='both'))
    ax3.set_xlabel('y')
    ax3.set_ylabel('z')
    ax3.set_aspect('equal', 'box')
    ax3.set_ylim(bottom=0)
    # cbar = fig.colorbar(cs2)


    cax = plt.axes([1.05, 0.11, 0.025, 0.78])
    cbar = fig.colorbar(
            cm.ScalarMappable(norm=kwargs['norm'], cmap=kwargs['cmap']),
            extend='both', cax = cax,
                            )

    # cbar.set_ticks(levels[::10])
    cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter(kwargs['tick_fmt']))
    cbar.set_label(kwargs['cbar_label'])

    fig.show()

    return fig, axes


def channel_two_plane(coords, data, **kwargs):
    import matplotlib.ticker as tick
    
    defaultKwargs = {
        'cmap'                  : cm.twilight_shifted,
        'figsize'               : (8, 6),
        'levels'                : 101,
        'norm'                  : colors.Normalize(data.min(), data.max()),
        'cbar_label'            : None,
        'dpi'                   : 150,
        'cbar_ticks'            : None,
        'yind'                  : None,
        'zind'                  : None,
    }


    kwargs = { **defaultKwargs, **kwargs}

    index = np.unravel_index(data.argmax(), data.shape)
    print(index)
    if kwargs['yind'] == None: 
        yind = index[1]
    else:
        yind = kwargs['yind']
        
    if kwargs['zind'] == None: 
        zind = index[2]
    else:
        zind = kwargs['zind']

    [x, y, z] = coords

    fig, axes = plt.subplots(4, figsize=kwargs['figsize'], dpi=kwargs['dpi'])

    gs = axes[1].get_gridspec()
    # remove the underlying axes
    for ax in axes[1:]:
        ax.remove()
    axbig = fig.add_subplot(gs[1:])

    print("x-z plane on y = %.4f" % y[yind])
    axes[0].contourf(x, z, data[:, yind, :].T, norm=kwargs['norm'], cmap=kwargs['cmap'], levels = kwargs['levels'])
    axes[0].set_aspect('equal', 'box')
    axes[0].set_ylim(0, 1)
    axes[0].set_xlim(0, 2*np.pi)
    axes[0].set_xticks([0,  np.pi, 2*np.pi], ["$0$", "$\pi$", "$2\pi$"])
    axes[0].set_ylabel("$z$")

    print("x-y plane on z = %.4f" % z[zind])
    axbig.contourf(x, y, data[:, :, zind].T, norm=kwargs['norm'], cmap=kwargs['cmap'], levels = kwargs['levels'])
    axbig.set_aspect('equal', 'box')
    axbig.set_ylim(0, np.pi)
    axbig.set_yticks([0, np.pi/2, np.pi], ["$0$", "$\pi/2$", "$\pi$"])
    axbig.set_xlim(0, 2*np.pi)
    axbig.set_xticks([0,  np.pi, 2*np.pi], ["$0$", "$\pi$", "$2\pi$"])
    axbig.set_xlabel("$x$")
    axbig.set_ylabel("$y$")

    fig.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.88, 0.11, 0.025, 0.78])
    cbar = fig.colorbar(
            cm.ScalarMappable(norm=kwargs['norm'], cmap=kwargs['cmap']),
            cax = cax,)
    if kwargs['cbar_ticks'] != None:
        cbar.set_ticks(kwargs['cbar_ticks'])
    cbar.formatter.set_powerlimits((0, 0))

    # to get 10^3 instead of 1e3
    cbar.formatter.set_useMathText(True)
    cbar.set_label(kwargs['cbar_label'])

    return fig, [axes[0], axbig, cax]


class plot_format():
    def __init__(self,):
        return
    
    
if __name__ == "__main__":
    1