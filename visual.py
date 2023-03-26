import numpy as np
import matplotlib.pyplot as plt
import scipy
import os
import sys
import glob
import contextlib
from PIL import Image
from generate_field import read_array_from_file, write_array_to_file

def write_array_to_file(filename, array,):
    out = []
    nx, ny, nz = array.shape

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                out.append(array[i, j, k])
    
    array_column_major = np.array(out)
    # Flatten the array in column-major order
    # array_column_major = np.reshape(array, (array.shape[2], array.shape[1], array.shape[0]), order='F').flatten()

    # Write the flattened array to a binary file
    with open(filename, "wb") as f:
        array_column_major.tofile(f)
    return array_column_major

def read_array_from_file(filename, shape, dtype=np.float64):
    # Read the flattened array from the binary file
    with open(filename, "rb") as f:
        array_column_major = np.fromfile(f, dtype=dtype)

    # Reshape the array into its original shape
    array = np.reshape(array_column_major, shape[::-1], order='C').T
    # array = np.reshape(array, shape[::-1], order='C')

    return array


def XYZ(domain, dims):
    # Define the range and number of points in each direction
    x_range = (0, domain[0])
    y_range = (0, domain[1])
    z_range = (0, domain[2])
    n_x, n_y, n_z = dims

    # Generate 1D arrays of x, y, and z coordinates
    x_coords = np.linspace(x_range[0], x_range[1], n_x)
    y_coords = np.linspace(y_range[0], y_range[1], n_y)
    z_coords = np.linspace(z_range[0], z_range[1], n_z)

    # Use meshgrid to generate a 3D grid of x, y, and z coordinates
    xx, yy, zz = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')

    # Display the shape of the resulting coordinate arrays
    print(f"Shape of xx: {xx.shape}")
    print(f"Shape of yy: {yy.shape}")
    print(f"Shape of zz: {zz.shape}")

    return xx, yy, zz

def png2gif(filename, varname):

    ROOT_DIR = os.path.abspath(filename)
    SLURM_ID = os.path.basename(filename)
    # filepaths
    fp_in = ROOT_DIR + "/%s.*.png" % varname
    fp_out = ROOT_DIR + "/%s.gif" % SLURM_ID

    # use exit stack to automatically close opened images
    with contextlib.ExitStack() as stack:

        # lazily load images
        imgs = (stack.enter_context(Image.open(f))
                for f in sorted(glob.glob(fp_in)))

        # extract  first image from iterator
        img = next(imgs)

        # https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif
        img.save(fp=fp_out, format='GIF', append_images=imgs,
                save_all=True, duration=200, loop=0)

class visual():
    """
    Visualize lesgo output
    """
    def __init__(self, dir) -> None:
        # Setup the Mesh
        self._setup(self, )

        # Lesgo outputs folder path
        self.dir = dir

    def _setup(self, domain = [2*np.pi, np.pi, 1], dims=[256, 256, 128], timesteps=1e4, gapT=500, out_format="%.4i",):
        self.dims = dims
        self.domain = domain
        self.nx, self.ny, self.nz = dims
        self.Lx, self.Ly, self.Lz = domain

        self.x_coord = np.linspace(0, self.Lx, self.nx)
        self.y_coord = np.linspace(0, self.Ly, self.ny)
        self.z_coord = np.linspace(0, self.Lz, self.nz)

        self.xx, self.yy, self.zz = XYZ(domain, dims)

        self.total_tau = timesteps
        self.dtau = gapT
        # Lesgo Output format of timestep
        self.oformat = out_format

        print('Mesh is created on domain %s with resolution %s' %(domain, dims))
        return 
    
    def plot_theta(self, path):
        t = 0
        theta_filename = self.dir+"/theta.{timestep:"+self.oformat+"}"
        while t<self.total_tau:
            theta = readfile(theta_filename.format(timestep = t), self.dims)   
        return
    

    def plume(self, data, ):
        # Create a 3D plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Define the X and Y arrays
        X, Z = np.meshgrid(self.x_coord, self.z_coord, )

        # Plot the surface
        ax.plot_surface(X, Z, data)

        # Add a contour to the plot
        ax.contour(X, Z, data, offset=-2, cmap='coolwarm')

        # Set the plot limits and labels
        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        ax.set_zlim(-2, 2)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        # Show the plot
        plt.show()
        return fig

    


if __name__ == "__main__":
    domain = np.array([2*np.pi, np.pi, 1])
    dims = [256, 256, 128]

    X, Y, Z = XYZ(domain, dims)




