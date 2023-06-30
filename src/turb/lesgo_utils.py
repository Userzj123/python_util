import numpy as np
import matplotlib.pyplot as plt
import pyutils.plot_utils as pltutils
import os

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
        defaultKwargs = {'domain': (2*np.pi, np.pi, 1), 'dims': (128, 128, 64) }
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

class lesgo_data():
    """
    Read lesgo's simulation result files.
    
    """
    
    def __init__(self, domain, dims, root_dir, ntheta=1) -> None:
        if not os.path.exists(root_dir):
            os.makedirs(root_dir)
            print('create root folder in %s' % root_dir)
        self.inputs_dir = root_dir + '/inputs'
        if not os.path.exists(self.inputs_dir):
            os.makedirs(self.inputs_dir)
            print('create inputs folder in %s' % self.inputs_dir)
            
        
        self.outputs_dir = root_dir + '/outputs'
        self.coords = coords_xyz(domain, dims, center=True, stretch=True)
        self.domain = domain
        self.dims = dims
        self.ntheta = ntheta
        
        self.data = {}
        
        self.ihalf_dims = dims*1
        self.ihalf_dims[0] += 1
        self.jhalf_dims = dims*1
        self.jhalf_dims[1] += 1
        
        self.ihalf_coords = coords_xyz(domain, self.ihalf_dims, center=False, stretch=True)
        self.jhalf_coords = coords_xyz(domain, self.jhalf_dims, center=False, stretch=True)
        
        # initialize dataset file names
        self.set_adjoint(adjoint=False)
        
    # Post-processing
    def read_data(self, t_ind):
        self.read_velocity(t_ind)
        self.read_scalar(t_ind)
        
        return
    
    def _fnames(self,):
        if not self.adjoint :
            self.u_f = self.outputs_dir + r'/baseflow/u_base.%.8i'
            self.v_f = self.outputs_dir + r'/baseflow/v_base.%.8i'
            self.w_f = self.outputs_dir + r'/baseflow/w_base.%.8i'
        else:
            self.u_f = self.outputs_dir + r'/u_velocity.%.8i'
            self.v_f = self.outputs_dir + r'/v_velocity.%.8i'
            self.w_f = self.outputs_dir + r'/w_velocity.%.8i'
        
        self.theta_f = self.outputs_dir + r'/theta.%.2i.%.8i'
    
    def set_adjoint(self, adjoint = False):
        self.adjoint = adjoint
        self._fnames()
        return
    
    def read_debug(self, t_ind):
        self.debug_advection_scalar(t_ind)
        self.debug_diffusion_scalar(t_ind)
        
        return
    
    def read_velocity(self, t_ind):
        u = read_array_from_file(self.u_f % t_ind, self.dims)
        v = read_array_from_file(self.v_f % t_ind, self.dims)
        w = read_array_from_file(self.w_f % t_ind, self.dims)
        
        self.data['u'] = u
        self.data['v'] = v
        self.data['w'] = w
        
        return u, v, w
    
    def _read_multi_scalar(self, data_f, t_ind, dims):
        data = read_array_from_file(data_f % (1, t_ind), dims)[np.newaxis, :, :, :]
        for k in range(1, self.ntheta):
            data = np.concatenate((data, read_array_from_file(data_f % (k+1, t_ind), dims)[np.newaxis, :, :, :]))
        return data
    
    
    def read_scalar(self, t_ind):
        thetas = self._read_multi_scalar(self.theta_f, t_ind, self.dims)
        
        self.data['theta'] = thetas
        return thetas
    
    def debug_advection_scalar(self, t_ind):
        adv_f = self.outputs_dir + r'/advection.%.2i.%.8i'
        
        dTdx_f = self.outputs_dir + r'/dTdx.%.2i.%.8i'
        dTdy_f = self.outputs_dir + r'/dTdy.%.2i.%.8i'
        dTdz_f = self.outputs_dir + r'/dTdz.%.2i.%.8i'
        
        self.data['adv'] = self._read_multi_scalar(adv_f, t_ind, self.dims)

        self.data['dTdx'] = self._read_multi_scalar(dTdx_f, t_ind, self.ihalf_dims)
        self.data['dTdy'] = self._read_multi_scalar(dTdy_f, t_ind, self.jhalf_dims)
        self.data['dTdz'] = self._read_multi_scalar(dTdz_f, t_ind, self.dims)

        return

    def debug_diffusion_scalar(self, t_ind):
        diff_f = self.outputs_dir + r'/diffusion.%.2i.%.8i'
        
        u_ihalf_f = self.outputs_dir + r'/u_ihalf.%.2i.%.8i'
        v_jhalf_f = self.outputs_dir + r'/v_jhalf.%.2i.%.8i'
        w_kw_f = self.outputs_dir + r'/w_kw.%.2i.%.8i'
        
        theta_ihalf_f = self.outputs_dir + r'/theta_ihalf.%.2i.%.8i'
        theta_jhalf_f = self.outputs_dir + r'/theta_jhalf.%.2i.%.8i'
        theta_kw_f = self.outputs_dir + r'/theta_kw.%.2i.%.8i'
        
        self.data["diff"] = self._read_multi_scalar(diff_f, t_ind, self.dims)
        
        self.data["u_ihalf"] = self._read_multi_scalar(u_ihalf_f, t_ind, self.ihalf_dims)
        self.data["v_jhalf"] = self._read_multi_scalar(v_jhalf_f, t_ind, self.jhalf_dims)
        self.data["w_kw"] = self._read_multi_scalar(w_kw_f, t_ind, self.dims)
        
        self.data["theta_ihalf"] = self._read_multi_scalar(theta_ihalf_f, t_ind, self.ihalf_dims)
        self.data["theta_jhalf"] = self._read_multi_scalar(theta_jhalf_f, t_ind, self.jhalf_dims)
        self.data["theta_kw"] = self._read_multi_scalar(theta_kw_f, t_ind, self.dims)
        
        return
    
    # [WorkingOn]
    # def sensor_location(self, loc_coords):
    #     """_summary_

    #     Args:
    #         loc_coords (_type_): input the location of sensor in Cartesian coordinate
            
    #     return: 
    #         coords_index (_type_): return the indices of sensor location in given domain.
    #     """
    #     self.npoints = len(loc_coords[0])
        
    #     # coords = [x, y, z], loc_coords = [px, py, pz]
    #     return
    
    # def sensor_measurements(self, tspan):
    #     tmin, tmax = tspan
    #     for t in range(int(tmin), int(tmax)):
    #         self.read_data(t)
    #         for ind, x in self.sensor_coords:
    #             measurement = self.data['theta'][:, x[0, ind], x[1, ind], x[2, ind]]
        
        
    #     return 
    
    
    # Pre-processing
        
    def generate_ic(self, **kwargs):
        defaultKwargs = {
            'ic_type' : 'gaussian',
            'homogeneous' : 'none',
        }
        kwargs = { **defaultKwargs, **kwargs }
        
        if kwargs['ic_type'] == 'gaussian':
            self.gaussian_ic(**kwargs)
            
    def gaussian_ic(self, **kwargs):
        """
        Generate a 3 dimensional numpy array with shape of [nx, ny, nz], which has a normal distribution with mean of mu_L=[x0, y0, z0] and variance of sigma. And the maximum value is set to be 1 in order to keep the magnitude same for differet IC.
        theta = e^(-((x-x0)^2 + (y-y0)^2 + (z-z0)^2)/(2*sigma^2))
        
        Parameters:
            domain : tuple of ints [Lx, Ly, Lz]
            dims   : tuple of ints [nx, ny, nz]
            mu_L   : tuple of ints [x0, y0, z0]
            sigma  : float, variance of the gaussian distribution
        Returns:
            data   : 3d numpy array with shape of [nx, ny, nz]
        
        """
        defaultKwargs = {
            'variance'          : 1e-1,
            'dtype'             : np.float64,
            'homogeneous'       : 'None',
            'source_point'      : [1, 1, 0.5],
            'ic_dir'            : self.inputs_dir,
            'varname'           : 'theta',
            'nk'                : 1,
            'readme_text'       : 'Domain =' + str(self.domain) + ', Dims = ' + str(self.dims) + ', nk = ' + str(kwargs['nk']) + ', Source at ' + str(kwargs['source_point'])
        }
        kwargs = { **defaultKwargs, **kwargs }
        
        x, y, z = self.coords
        x0, y0, z0 = kwargs['source_point']
        
        x_homo = 0 if 'x' in kwargs['homogeneous'].lower() else 1
        y_homo = 0 if 'y' in kwargs['homogeneous'].lower() else 1
        z_homo = 0 if 'z' in kwargs['homogeneous'].lower() else 1

        
        data = np.zeros(self.dims, dtype=kwargs['dtype'])
        for i, xx in enumerate(x):
            for j, yy in enumerate(y):
                for k, zz in enumerate(z):
                        data[i, j, k] = np.exp(-( (xx-x0)**2*x_homo + (yy-y0)**2*y_homo + (zz-z0)**2*z_homo )/kwargs['variance']**2/2)
        data = data/data.max()
        
        
        ic_fname = kwargs['ic_dir'] + '/%s.IC.%.2i' %(kwargs['varname'], kwargs['nk'])
        
        ic_readme_fname = self.inputs_dir + '/readme.md'
        from datetime import datetime
            
        with open(ic_readme_fname, "a+") as f:
            f.write(datetime.today().strftime('%Y-%m-%d'))
            f.write('  Update Initial Condition\n')
            f.write(kwargs['readme_text'] + '\n')
        
        return write_array_to_file(ic_fname, data, **kwargs)
    
    def source_sameasic(self,):
        import shutil            
        ic_files = [filename for filename in os.listdir(self.inputs_dir) if filename.startswith("theta")]
        new_fname = ['./source.'+fname.split('.')[-1] for fname in ic_files]
        [shutil.copyfile(self.inputs_dir+ '/' + ic_files[ind], self.inputs_dir+name) for ind, name in enumerate(new_fname)]
        return new_fname
        
        
if __name__ == "__main__":
    result_dir = '/home/zyou6474/tasks/channel_flow'
    dims = (128, 128, 64)
    domain = (2*np.pi, np.pi, 1)

    ldata = lesgo_data(domain, dims, result_dir, ntheta=3)


    ldata.read_data(1)
    ldata.data['theta'].shape
