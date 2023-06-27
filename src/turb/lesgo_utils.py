import numpy as np
import matplotlib.pyplot as plt
import pyutils.plot_utils as pltutils

def write_array_to_file(filename, array, gen_fig=False, **kwargs):
    """
    Write input numpy arrray with dims [nx, ny, nz], which is written and read in row-major order (first index - nz changes fastest), into a binary file.
    Parameters:
        filename : String
        array    : 3d numpy array with shape [nx, ny, nz]
    Returns:
        array_column_major : 1d numpy array with shape [nx*ny*nz], written in column-major order(last index - nx changes fastest)
    
    """
    array_column_major = np.reshape(array, (-1), order='F')

    with open(filename, "wb") as f:
        array_column_major.tofile(f)

    if gen_fig:
        defaultKwargs = {'domain': (2*np.pi, 2*np.pi, 1), 'dims': (128, 128, 64) }
        kwargs = { **defaultKwargs, **kwargs }
        
        coords = coords_xyz(kwargs['domain'], kwargs['dims'])
        fig = pltutils.contour_channel(coords, array)
        return fig, array_column_major
    else:
        return array_column_major
    
def read_array_from_file(filename, shape, dtype=np.float64):
    """
    Read the flattened array from the binary file
    Parameters:
        filename : String
        shape    : int or tuple of ints [nx, ny, nz]
    Returns:
        array    : Numpy array with shape [nx, ny, nz] read from binary file.
    """
    with open(filename, "rb") as f:
        array_column_major = np.fromfile(f, dtype=dtype)

    # Reshape the array into its original shape
    array = np.reshape(array_column_major, shape[::-1], order='C').T
    # array = np.reshape(array, shape[::-1], order='C')

    return array

def coords_xyz(domain, dims, stretch=False, str_factor=1.5, center=False):
    # Define the range and number of points in each direction
    x_range = (0, domain[0])
    y_range = (0, domain[1])
    z_range = (0, domain[2])
    n_x, n_y, n_z = dims

    if center:
        # Generate 1D arrays of x, y, and z coordinates
        x_coords = np.linspace(x_range[0], x_range[1], int(n_x)+1)
        y_coords = np.linspace(y_range[0], y_range[1], int(n_y)+1)
        z_coords = np.linspace(z_range[0], z_range[1], int(n_z)+1)
    
        x_coords = 0.5 * x_coords[:-1] + 0.5 * x_coords[1:]
        y_coords = 0.5 * y_coords[:-1] + 0.5 * y_coords[1:]
        z_coords = 0.5 * z_coords[:-1] + 0.5 * z_coords[1:]
        
    else:
        # Generate 1D arrays of x, y, and z coordinates
        x_coords = np.linspace(x_range[0], x_range[1], int(n_x))
        y_coords = np.linspace(y_range[0], y_range[1], int(n_y))
        z_coords = np.linspace(z_range[0], z_range[1], int(n_z))

    # k(uv) grid
    if stretch:
        z_stretch = z_range[1]*(1+(np.tanh(str_factor*(z_coords/z_range[1]-1))/np.tanh(str_factor)))
        return x_coords, y_coords, z_stretch
    
    return x_coords, y_coords, z_coords

class lesgo():
    """
    
    """
    
    def __init__(self, domain, dims, result_dir) -> None:
        self.dir = result_dir
        self.coords = coords_xyz(domain, dims)
        self.domain = domain
        self.dims = dims
        
    
    def read_velocity(self, t_ind):
        u_f = self.dir + r'/u.%.8i'
        v_f = self.dir + r'/v.%.8i'
        w_f = self.dir + r'/w.%.8i'
        
        u = read_array_from_file(u_f % t_ind, self.dims)
        v = read_array_from_file(v_f % t_ind, self.dims)
        w = read_array_from_file(w_f % t_ind, self.dims)

        return u, v, w
    
    def read_scalar(self, t_ind, nk):
        theta_f = self.dir + r'/theta.%.2i.%.8i'
        
        theta = read_array_from_file(theta_f % (nk, t_ind), self.dims)
        return theta
    
    def debug_advection_scalar(self, t_ind, nk):
        adv_f = self.dir + r'/advection.%.2i.%.8i'
        
        dTdx_f = self.dir + r'/dTdx.%.2i.%.8i'
        dTdy_f = self.dir + r'/dTdy.%.2i.%.8i'
        dTdz_f = self.dir + r'/dTdz.%.2i.%.8i'
        
        data = {}
        data['adv'] = read_array_from_file(adv_f % (nk, t_ind), self.dims)

        return data

    def debug_diffusion_scalar(self, t_ind, nk):
        diff_f = self.dir + r'/diffusion.%.2i.%.8i'
        
        u_ihalf_f = self.dir + r'/u_ihalf.%.2i.%.8i'
        v_jhalf_f = self.dir + r'/v_jhalf.%.2i.%.8i'
        w_kw_f = self.dir + r'/w_kw.%.2i.%.8i'
        
        theta_ihalf_f = self.dir + r'/theta_ihalf.%.2i.%.8i'
        theta_jhalf_f = self.dir + r'/theta_jhalf.%.2i.%.8i'
        theta_kw_f = self.dir + r'/theta_kw.%.2i.%.8i'
        
        data = {}
        data["diff"] = read_array_from_file(diff_f % (nk, t_ind), self.dims)
        
        data["u_ihalf"] = read_array_from_file(u_ihalf_f % (nk, t_ind), self.dims)
        data["v_jhalf"] = read_array_from_file(v_jhalf_f % (nk, t_ind), self.dims)
        data["w_kw"] = read_array_from_file(w_kw_f % (nk, t_ind), self.dims)
        
        data["theta_ihalf"] = read_array_from_file(theta_ihalf_f % (nk, t_ind), self.dims)
        data["theta_jhalf"] = read_array_from_file(theta_jhalf_f % (nk, t_ind), self.dims)
        data["theta_kw"] = read_array_from_file(theta_kw_f % (nk, t_ind), self.dims)
        
        return data
        
