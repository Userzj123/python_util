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

def plot_ins_3d(path, dims):

    with open(path, 'rb') as f:
        U = np.fromfile(f, dtype=np.float64,)

    U = np.reshape(U, dims)

    # Define dimensions
    Nz, Ny, Nx = dims
    Z, Y, X = np.meshgrid(np.arange(Nz), np.arange(Ny), np.arange(Nx), indexing="ij")

    kw = {
        'vmin': U.min(),
        'vmax': U.max(),
        'levels': np.linspace(U.min(), U.max(), 10),
    }

    # Create a figure with 3D ax
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111, projection='3d')

    # Plot contour surfaces
    _ = ax.contourf(
        X[:, :, 0], Y[:, :, 0], U[:, :, 0],
        zdir='z', offset=0, **kw
    )
    _ = ax.contourf(
        X[0, :, :], U[0, :, :], Z[0, :, :],
        zdir='y', offset=0, **kw
    )
    C = ax.contourf(
        U[:, -1, :], Y[:, -1, :], Z[:, -1, :],
        zdir='x', offset=X.max(), **kw
    )
    # --


    # Set limits of the plot from coord limits
    xmin, xmax = X.min(), X.max()
    ymin, ymax = Y.min(), Y.max()
    zmin, zmax = Z.min(), Z.max()
    ax.set(xlim=[xmin, xmax], ylim=[ymin, ymax], zlim=[zmin, zmax])

    # Plot edges
    edges_kw = dict(color='0.4', linewidth=1, zorder=1e3)
    ax.plot([xmax, xmax], [ymin, ymax], 0, **edges_kw)
    ax.plot([xmin, xmax], [ymin, ymin], 0, **edges_kw)
    ax.plot([xmax, xmax], [ymin, ymin], [zmin, zmax], **edges_kw)

    # Set labels and zticks
    ax.set(
        xlabel='X [km]',
        ylabel='Y [km]',
        zlabel='Z [m]',
        zticks=[0, -150, -300, -450],
    )

    # Set zoom and angle view
    ax.view_init(40, -30, 0)
    ax.set_box_aspect(None, zoom=0.9)

    # Colorbar
    fig.colorbar(C, ax=ax, fraction=0.02, pad=0.1, label='Name [units]')

    # Show Figure
    fig.savefig(r'/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/test.png')
    return 1


def plot_instantaneous(path, dims):
    U = read_array_from_file(path, dims)

    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    ax.contourf(U[:, :, 64])
    ax.set_aspect('equal', 'box')
    # ax[1].contourf(U)
    # ax[1].set_aspect('equal', 'box')
    fig.savefig(r'/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/test.png')
    return fig


def contourf_t(f, t, domain, dims, z=64):

    theta = read_array_from_file(f % t, dims)
    # theta[theta<theta.max()*0.01]=0
    levs = 10**np.linspace(-6, 0, 101)

    # theta[theta<0.1*theta.max()] = 0
    
    # Generate the coordinates
    x_coords, y_coords, z_coords = xyz(domain, dims)

    fig, ax = plt.subplots(figsize=(8,4), dpi=150)
    cs = ax.contourf(x_coords, y_coords, theta[:, :, z].T, levs, norm=colors.LogNorm(), cmap=cm.binary)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.axis('equal')

    ax.set_ylim(bottom=0)

    cbar = fig.colorbar(cs)
    cbar.set_ticks(levs[::50])
    cbar.set_label(r'Scalar Concentration', size=16)

    ax.set_title('z=0.5')
    # fig.show()
    return fig

def output_theta():
    dir = '/scratch4/qwang4/ext-zyou6474/channel_flow_multi_scalar.14107081'
    domain = [2*np.pi, np.pi, 1]
    dims = [256, 256, 128]

    theta_f = dir + '/outputs/theta.01.%.8i'
    u_f = dir + '/outputs/baseflow/u_base.%.8i'
    v_f = dir + '/outputs/baseflow/v_base.%.8i'
    w_f = dir + '/outputs/baseflow/w_base.%.8i'
    
    filenames = []
    timestep = np.linspace(0, 10000, 101)
    for i in timestep:
        fig = contourf_t(theta_f, i, domain, dims)

        # create file name and append it to a list
        filename = f'results/%.5i.png' % i
        filenames.append(filename)

        # save frame
        fig.savefig(filename)
        plt.close()
    # build gif
    with imageio.get_writer('results/mygif.gif', mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Remove files
    for filename in set(filenames):
        os.remove(filename)


if __name__ == "__main__":
    output_theta()
