import numpy as np 
import matplotlib.pyplot as plt


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
    with open(path, 'rb') as f:
        U = np.fromfile(f, dtype=np.float64,)

    U = np.reshape(U, dims)

    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    ax.contourf(U[20, :, :])
    ax.set_aspect('equal', 'box')
    # ax[1].contourf(U)
    # ax[1].set_aspect('equal', 'box')
    fig.savefig(r'/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/test.png')
    return fig





if __name__ == "__main__":
    dims = [64, 128, 128]
    dims = np.array(dims)*2
    plot_instantaneous('/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/inputs/theta.00000000', dims,)