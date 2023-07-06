import numpy as np
import matplotlib.pyplot as plt
import pyutils.plot_utils as pltutils
from pyutils.cartesian import coords_xyz, meshgrid
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
    print('write data into %s' % filename)

    if gen_fig:
        defaultKwargs = {'domain': (2*np.pi, np.pi, 1), 'dims': (128, 128, 64) }
        kwargs = { **defaultKwargs, **kwargs }
        
        coords = coords_xyz(kwargs['domain'], kwargs['dims'], center=True, stretch=True)
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
        self.ntheta : np.int64 = ntheta
        
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
    
    def _fnames(self, **fmt_kwargs):
        defaultKwargs = {
            'fmt_tstep'   : r'%.8i',
            'fmt_ntheta'   : r'%.2i'
        }
        self.fmt_kwargs = { **defaultKwargs, **fmt_kwargs }
        
        if not self.adjoint :
            self.u_f = self.outputs_dir + r'/baseflow/u_base.' + self.fmt_kwargs['fmt_tstep']
            self.v_f = self.outputs_dir + r'/baseflow/v_base.' + self.fmt_kwargs['fmt_tstep']
            self.w_f = self.outputs_dir + r'/baseflow/w_base.' + self.fmt_kwargs['fmt_tstep']
        else:
            self.u_f = self.outputs_dir + r'/u_velocity.' + self.fmt_kwargs['fmt_tstep']
            self.v_f = self.outputs_dir + r'/v_velocity.' + self.fmt_kwargs['fmt_tstep']
            self.w_f = self.outputs_dir + r'/w_velocity.' + self.fmt_kwargs['fmt_tstep']
        
        self.theta_f = self.outputs_dir + r'/theta.' + self.fmt_kwargs['fmt_ntheta'] + '.' + self.fmt_kwargs['fmt_tstep']
    
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
        adv_f = self.outputs_dir + r'/advection.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        
        dTdx_f = self.outputs_dir + r'/dTdx.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        dTdy_f = self.outputs_dir + r'/dTdy.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        dTdz_f = self.outputs_dir + r'/dTdz.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        
        self.data['adv'] = self._read_multi_scalar(adv_f, t_ind, self.dims)

        self.data['dTdx'] = self._read_multi_scalar(dTdx_f, t_ind, self.ihalf_dims)
        self.data['dTdy'] = self._read_multi_scalar(dTdy_f, t_ind, self.jhalf_dims)
        self.data['dTdz'] = self._read_multi_scalar(dTdz_f, t_ind, self.dims)

        return

    def debug_diffusion_scalar(self, t_ind):
        diff_f = self.outputs_dir + r'/diffusion.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        
        u_ihalf_f = self.outputs_dir + r'/u_ihalf.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        v_jhalf_f = self.outputs_dir + r'/v_jhalf.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        w_kw_f = self.outputs_dir + r'/w_kw.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        
        theta_ihalf_f = self.outputs_dir + r'/theta_ihalf.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        theta_jhalf_f = self.outputs_dir + r'/theta_jhalf.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        theta_kw_f = self.outputs_dir + r'/theta_kw.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        
        self.data["diff"] = self._read_multi_scalar(diff_f, t_ind, self.dims)
        
        self.data["u_ihalf"] = self._read_multi_scalar(u_ihalf_f, t_ind, self.ihalf_dims)
        self.data["v_jhalf"] = self._read_multi_scalar(v_jhalf_f, t_ind, self.jhalf_dims)
        self.data["w_kw"] = self._read_multi_scalar(w_kw_f, t_ind, self.dims)
        
        self.data["theta_ihalf"] = self._read_multi_scalar(theta_ihalf_f, t_ind, self.ihalf_dims)
        self.data["theta_jhalf"] = self._read_multi_scalar(theta_jhalf_f, t_ind, self.jhalf_dims)
        self.data["theta_kw"] = self._read_multi_scalar(theta_kw_f, t_ind, self.dims)
        
        return
    
    # [Done] - Jul 5, 2023
    def sensor_init_(self, s_locs):
        """_summary_

        Args:
            loc_coords (_type_): input the location of sensor in Cartesian coordinate
            
        return: 
            coords_index (_type_): return the indices of sensor location in given domain.
        """
        self.nsensors = len(s_locs[0])
        
        # coords = [x, y, z], loc_coords = [px, py, pz]
        
        self.s_locs = s_locs
        return
    
    def channel_obs(self, **kwargs):
        from matplotlib import cm, colors, ticker
        import matplotlib.ticker as tick
        from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
        defaultKwargs = {
            'alpha' : 0.1,
            'vmin' : self.data['theta'].min(),
            'vmax' : self.data['theta'].max(),
            'cmap' : cm.Oranges,
            'levels': 101,
            'figsize' : (8, 8),
            'dpi'   : 150,
            'nk'    : 1,
            'y_ind' : 64,
            'clip'  : False,
            'norm'  : colors.Normalize
        }
        kwargs = { **defaultKwargs, **kwargs }
        
        levels = np.linspace(kwargs['vmin'], kwargs['vmax'], kwargs['levels'])

        fig, axes = plt.subplots(self.ntheta+1, 1, figsize=kwargs['figsize'], dpi = kwargs['dpi'])
        axes = axes.flatten()

        ct1 = axes[0].contourf(self.coords[0], self.coords[2], self.data['theta'][kwargs['nk'], :, kwargs['y_ind'], :].T, levels=kwargs['levels'],
                    norm=kwargs['norm'](vmin=kwargs['vmin'], vmax=kwargs['vmax'], clip=kwargs['clip']), cmap=kwargs['cmap'], 
                    )
        axes[0].scatter(self.s_locs[0], self.s_locs[-1], marker='o', color='black', s = 1)
        axes[0].set_xlabel('x')
        axes[0].set_ylabel('z')
        axes[0].set_aspect('equal', 'box')
        axes[0].set_ylim(bottom=0)
        
        ax1_divider = make_axes_locatable(axes[0])
        # Add an Axes to the right of the main Axes.
        cax1 = ax1_divider.append_axes("right", size="2%", pad="1%")
        cbar = fig.colorbar(
            cm.ScalarMappable(norm=kwargs['norm'](vmin=kwargs['vmin'], vmax=kwargs['vmax'], clip=kwargs['clip']), cmap=kwargs['cmap']),
            extend='both', cax = cax1,
                            )

        # cbar.set_ticks(levels[::10])
        # cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))

        for i in range(self.nsensors):
            axes[i+1].plot(self.obs_t, self.data['obs'][i, kwargs['nk'], :], color='black', label='obs', linewidth=0.5)
            axes[i+1].set_xlabel('t')
            axes[i+1].set_ylabel('theta')
            axes[i+1].set_ylim(kwargs['vmin'], kwargs['vmax'])
            axes[i+1].set_xlim(self.obs_t.min(), self.obs_t.max())
            axes[i+1].legend(frameon=False)
            
        return fig, axes
    
    
    def sensor_measurements(self, sensor_locs, tt, **kwargs):
        import imageio
        from datetime import datetime
        defaultKwargs = {
            'gif_fname' : './result' 
        }
        kwargs = { **defaultKwargs, **kwargs }
        
        kwargs['gif_fname'] += '_%s.gif' % datetime.today().strftime('%Y-%m-%d')
        
        self.sensor_init_(sensor_locs)
        self.obs_t = tt
        
        mesh = meshgrid(domain=self.domain, dims = self.dims)
        xi, yi, zi = mesh.cartesian2index(sensor_locs)
        
        gif_dir = './tmp_gif'
        if not os.path.exists(gif_dir):
            os.makedirs(gif_dir)
        filenames = []
    
    
        nsensor = len(xi)
        self.data['obs'] = np.full(shape=(nsensor, self.ntheta, len(tt)), fill_value=np.nan)
        for t_ind, t in enumerate(tt):
            self.read_scalar(t) 
            for ns in range(self.ntheta):
                for ind in range(nsensor):
                    self.data['obs'][ind, ns, t_ind] = self.data['theta'][ns, xi[ind], yi[ind], zi[ind]] 
                
            fig, ax = self.channel_obs(**kwargs)
            # create file name and append it to a list
            filename = gif_dir + f'/%.5i.png' % t
            filenames.append(filename)
            
            # save frame
            fig.savefig(filename, bbox_inches='tight')
            plt.close()
            
        # build gif
        with imageio.get_writer(kwargs['gif_fname'], mode='I') as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)
        print('gif written in %s' % kwargs['gif_fname'])
                
        # Remove files
        for filename in set(filenames):
            os.remove(filename)
            
        return
        
    def sensor_gaussian_ic(self, nsensors=3):
        
        return
    
    
    # Pre-processing
        
    def generate_ic(self, **kwargs):
        defaultKwargs = {
            'ic_type' : 'gaussian',
            'homogeneous' : 'none',
        }
        kwargs = { **defaultKwargs, **kwargs }
        
        if kwargs['ic_type'] == 'gaussian':
            self.gaussian_ic(**kwargs)
            
    def gaussian_field(self, **kwargs):
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
            'input_dir'         : self.inputs_dir,
            'fieldname'         : 'theta.IC',
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
        
        
        fname = (kwargs['input_dir'] + '/%s.' + self.fmt_kwargs['fmt_ntheta']) %(kwargs['fieldname'], kwargs['nk'])
        
        ic_readme_fname = self.inputs_dir + '/readme.md'
        from datetime import datetime
            
        with open(ic_readme_fname, "a+") as f:
            f.write(datetime.today().strftime('%Y-%m-%d'))
            f.write('Update %s\n' % fname)
            f.write(kwargs['readme_text'] + '\n')
        
        return write_array_to_file(fname, data, **kwargs)
    
              
    def constant_field(self, **kwargs):
        defaultKwargs = {
            'const'             : 0,
            'dtype'             : np.float64,
            'input_dir'         : self.inputs_dir,
            'fieldname'         : 'theta.IC',
            'nk'                : 1,
            'readme_text'        : ''
        }
        kwargs = { **defaultKwargs, **kwargs }
        kwargs['gen_fig'] = False
        
        kwargs['readme_text'] += 'Domain =' + str(self.domain) + ', Dims = ' + str(self.dims) + ', nk = ' + str(kwargs['nk']) + ', Constant of ' + str(kwargs['const'])


        
        data = np.zeros(self.dims, dtype=kwargs['dtype'])
        data += kwargs['const']
        
        
        fname = (kwargs['input_dir'] + '/%s.' + self.fmt_kwargs['fmt_ntheta']) %(kwargs['fieldname'], kwargs['nk'])
        
        ic_readme_fname = self.inputs_dir + '/readme.md'
        from datetime import datetime
            
        with open(ic_readme_fname, "a+") as f:
            f.write(datetime.today().strftime('%Y-%m-%d'))
            f.write('Update %s\n' % fname)
            f.write(kwargs['readme_text'] + '\n')
        
        return write_array_to_file(fname, data, **kwargs)
    
    def source_sameasic(self,):
        import shutil            
        ic_files = [filename for filename in os.listdir(self.inputs_dir) if filename.startswith("theta")]
        new_fname = ['/source.'+fname.split('.')[-1] for fname in ic_files]
        [shutil.copyfile(self.inputs_dir+ '/' + ic_files[ind], self.inputs_dir+name) for ind, name in enumerate(new_fname)]
        return new_fname
        
    # [Working On]
    def adjoint_ratio(self,):
        return
        
        
if __name__ == "__main__":
    root_dir = '/home/zyou6474/tasks/adjoint_steady_channel_flow'
    dims = [128, 128, 64]
    domain = [2*np.pi, np.pi, 1]

    ldata = lesgo_data(domain, dims, root_dir, ntheta=3)

    ldata.gaussian_field(fieldname = 'theta.IC', nk=3, source_point = [np.pi*5/3, np.pi/2, 2/3], gen_fig=True, variance=1e-2)