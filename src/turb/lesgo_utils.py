import numpy as np
import matplotlib.pyplot as plt
import pyutils.plot_utils as pltutils
from pyutils.cartesian import coords_xyz, meshgrid
import os
from configparser import ConfigParser

from datetime import datetime
date_marker = datetime.today().strftime('%Y_%m_%d')


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
    
    if 'readme_text' in kwargs.keys():
        ic_readme_fname = kwargs['input_dir'] + '/readme.md'
        
        with open(ic_readme_fname, "a+") as f:
            f.write(date_marker)
            f.write('Update %s\n' % filename)
            f.write(kwargs['readme_text'] + '\n')

    if gen_fig:
        defaultKwargs = {'domain': (2*np.pi, np.pi, 1), 'dims': (128, 128, 64) }
        kwargs = { **defaultKwargs, **kwargs }
        
        coords = coords_xyz(kwargs['domain'], kwargs['dims'], center=True, stretch=True)
        fig = pltutils.contour_channel(coords, array)
        return fig, array
    else:
        return array
    
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

def read_lesgoconf(lesgoconf_fname):        
    with open(lesgoconf_fname, 'r') as f:
        cmd = 0
        configs = {}
        block_tmp = {}
        block_entry = 0
        
        for line in f.readlines():
            if '!' in line:
                if cmd == 0:
                    comment = []
                    cmd  += 1
                    
                comment.append(line[line.index("!"):])
                line = line[:line.index("!")]
            else:
                cmd = 0
            
            if line.strip():
                if block_entry == 1 and '}' not in line:
                    varname, value = line.split('=')
                    varname = varname.strip()
                    values = value.replace('\n', '').strip().split()
                    
                    # if len(values) == 1:
                    #     if values[0] == '.true.':
                    #         block_tmp[varname] = True
                    #     elif values[0] == '.false.':
                    #         block_tmp[varname] = False
                    #     else:
                    #         try: block_tmp[varname] = float(values[0])
                    #         except:  block_tmp[varname] = value
                    # else:
                    #     try: block_tmp[varname] = [float(d) for d in values]
                    #     except:  block_tmp[varname] = value
                    block_tmp[varname] = value.replace('\n', '')
                    
                if '{' in line:
                    contents = line.split('{')
                    block_name = contents[0].strip()
                    block_entry = 1
                elif '}' in line:
                    configs[block_name] = block_tmp
                    block_entry = 2
                    block_tmp = {}
    
    return configs

def write_lesgoconf(fname, config: ConfigParser):
    with open(fname, 'w') as f:
        for sec in config.sections():
            f.write(sec + '{ \n')
            for var in config[sec]:
                f.write(var + ' = ' + config[sec][var] + '\n')
                
            f.write('} \n\n')
    print('Write configs in %s' % fname)
    return 


def unit_gaussian(mesh: meshgrid, x0, var, homo_axes=''):
    """    Generate a 3 dimensional numpy array with shape of [nx, ny, nz], which has a normal distribution with mean of mu_L=[x0, y0, z0] and variance of sigma. And the maximum value is set to be 1 in order to keep the magnitude same for differet IC.
    theta = e^(-((x-x0)^2 + (y-y0)^2 + (z-z0)^2)/(2*sigma^2))

    Args:
        mesh (meshgrid): _description_
        x0 (_type_): _description_
        var (_type_): _description_
        homo_axes (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    
    x, y, z = mesh.coords
    x0, y0, z0 = x0
    
    x_homo = 0 if 'x' in homo_axes.lower() else 1
    y_homo = 0 if 'y' in homo_axes.lower() else 1
    z_homo = 0 if 'z' in homo_axes.lower() else 1

    
    data = np.zeros(mesh.shape)
    for i, xx in enumerate(x):
        for j, yy in enumerate(y):
            for k, zz in enumerate(z):
                    data[i, j, k] = np.exp(-( (xx-x0)**2*x_homo + (yy-y0)**2*y_homo + (zz-z0)**2*z_homo )/var**2/2)
    data = data/data.max()
    
    return data

class lesgo_data():
    """
    Read lesgo's simulation result files.
    
    """
    
    def __init__(self, domain, dims, forward_dir, ntheta=1, center=True, stretch=True) -> None:
        self.forward_dir = forward_dir
        if not os.path.exists(forward_dir):
            os.makedirs(forward_dir)
            print('create root folder in %s' % forward_dir)
        self.forward_inputs_dir = forward_dir + '/inputs'
        self.inputs_dir = self.forward_inputs_dir
        if not os.path.exists(self.inputs_dir):
            os.makedirs(self.inputs_dir)
            print('create inputs folder in %s' % self.inputs_dir)
            
        
        self.forward_output_dir = forward_dir + '/output'
        self.output_dir = self.forward_output_dir
        
        self.coords = coords_xyz(domain, dims, center=center, stretch=stretch)
        self.domain = domain
        self.dims = dims
        self.ntheta : np.int64 = ntheta
        
        self.data = {}
        
        self.ihalf_dims = dims*1
        self.ihalf_dims[0] += 1
        self.jhalf_dims = dims*1
        self.jhalf_dims[1] += 1
        self.kw_dims = dims*1
        self.kw_dims[2] += 1
        
        self.ihalf_coords = coords_xyz(domain, self.ihalf_dims, center=False, stretch=True)
        self.jhalf_coords = coords_xyz(domain, self.jhalf_dims, center=False, stretch=True)
        self.kw_coords = coords_xyz(domain, self.kw_dims, center=False, stretch=True)
        
        # initialize dataset file names
        self.set_adjoint(adjoint=False)
        
        self.source_coords = np.empty(shape=(0, 3))
        
    # Post-processing
    
    
    def read_data(self, t_ind):
        self.read_velocity(t_ind)
        self.read_scalar(t_ind)
        
        return
    
    def read_velocity(self, t_ind, f=None):
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
    
    def read_scalar(self, t_ind, threshold=None):
        thetas = self._read_multi_scalar(self.theta_f, t_ind, self.dims)
        
        if threshold != None:
            thetas[thetas<threshold] = 0
        
        self.data['theta'] = thetas
        return thetas
    
    def read_inputs(self, dir):
        """
        Read the data of inputs field.
        """
        if not self.adjoint :
            u_icf = dir + '/u_velocity.IC'
            v_icf = dir + '/v_velocity.IC'
            w_icf = dir + '/w_velocity.IC'
            
            self.data['u_ic'] = read_array_from_file(u_icf, self.dims)
            self.data['v_ic'] = read_array_from_file(v_icf, self.dims)
            self.data['w_ic'] = read_array_from_file(w_icf, self.dims)
        
        source_fs = [os.path.join(dir, filename) for filename in os.listdir(dir) if filename.startswith("source.")]
        theta_IC_fs = [os.path.join(dir, filename) for filename in os.listdir(dir) if filename.startswith("theta.")]

        self.data['source'] = np.empty((0, self.dims[0], self.dims[1], self.dims[2]))
        for find, fn in enumerate(source_fs):
            self.data['source'] = np.concatenate((self.data['source'], read_array_from_file(fn, self.dims)[np.newaxis, :, :, :]))
        
        
        self.data['theta_IC'] = np.empty((0, self.dims[0], self.dims[1], self.dims[2]))
        for find, fn in enumerate(theta_IC_fs):
            self.data['theta_IC'] = np.concatenate((self.data['theta_IC'], read_array_from_file(fn, self.dims)[np.newaxis, :, :, :]))
        
        return
    
    def write_inputs(self, ):
        # Velocity Initial Conditions
        write_array_to_file(self.inputs_dir + '/u_velocity.IC', self.data['u_ic'])
        write_array_to_file(self.inputs_dir + '/v_velocity.IC', self.data['v_ic'])
        write_array_to_file(self.inputs_dir + '/w_velocity.IC', self.data['w_ic'])
        
        # Scalar Initial Conditions
        self.ntheta = self.data['theta_IC'].shape[0]
        self.config['SCALAR']['ntheta'] = (self.fmt_kwargs['fmt_ntheta']) % self.ntheta
        thetaic_f = self.inputs_dir + '/theta.' + self.fmt_kwargs['fmt_ntheta'] + '.IC'
        for nk in range(self.ntheta):
            write_array_to_file(thetaic_f % (nk+1), self.data['theta_IC'][nk])
            
        # Source
        if self.config['SCALAR']['source_opt'] == '1':
            if self.data['source'].shape[0] != self.ntheta:
                raise Exception("Sources of all scalar need to be defined in self.data['source'] with shape of (nk, nx, ny, nz)")
            else:   
                source_f = self.inputs_dir + '/source.' + self.fmt_kwargs['fmt_ntheta']
                for nk in range(self.ntheta):
                    write_array_to_file(source_f % (nk+1), self.data['source'][nk])
            
        # LESGO.CONF
        ## Update shape of domain 
        self.config['DOMAIN']['nx'] = f'{self.dims[0]}'
        self.config['DOMAIN']['ny'] = f'{self.dims[1]}'
        self.config['DOMAIN']['nz'] = f'{self.dims[2]}'
        self.config['DOMAIN']['lx'] = f'{self.domain[0]}'
        self.config['DOMAIN']['ly'] = f'{self.domain[1]}'
        self.config['DOMAIN']['lz'] = f'{self.domain[2]}'
        ## Write lesgo.conf
        write_lesgoconf(self.inputs_dir + '/lesgo.conf', self.config)
        
        return 

    def _fnames(self, **fmt_kwargs):
        defaultKwargs = {
            'fmt_tstep'   : r'%.8i',
            'fmt_ntheta'   : r'%.2i'
        }
        self.fmt_kwargs = { **defaultKwargs, **fmt_kwargs }
        

        self.u_f = self.output_dir + r'/velocity/u_velocity.' + self.fmt_kwargs['fmt_tstep']
        self.v_f = self.output_dir + r'/velocity/v_velocity.' + self.fmt_kwargs['fmt_tstep']
        self.w_f = self.output_dir + r'/velocity/w_velocity.' + self.fmt_kwargs['fmt_tstep']

        
        self.theta_f = self.output_dir + r'/scalar/theta.' + self.fmt_kwargs['fmt_ntheta'] + '.' + self.fmt_kwargs['fmt_tstep']
    
        if self.adjoint:
            self.adjoint_f = self.adjoint_dir + r'/scalar/theta.' + self.fmt_kwargs['fmt_ntheta'] + '.' + self.fmt_kwargs['fmt_tstep']
            
        return 
    
    
    def set_adjoint(self, adjoint = False, adjoint_dir = None, **kwargs):
        self.adjoint = adjoint
        if adjoint:
            self.adjoint_dir = adjoint_dir
            self.adjoint_inputs_dir = self.adjoint_dir + '/inputs'
            self.adjoint_output_dir = self.adjoint_dir + '/output'
            
            self.inputs_dir = self.adjoint_inputs_dir
            self.output_dir = self.adjoint_output_dir
            if not os.path.exists(adjoint_dir):
                os.makedirs(adjoint_dir)
                print('create adjoint root folder in %s' % adjoint_dir)
            if not os.path.exists(self.adjoint_inputs_dir):
                os.makedirs(self.adjoint_inputs_dir)
                print('create adjoint inputs folder in %s' % self.adjoint_inputs_dir)
        else:
            self.inputs_dir = self.forward_inputs_dir
            self.output_dir = self.forward_output_dir
                
        self._fnames(**kwargs)
        return
    
    # [Working on] maybe combine adjoint and forward init together as well as previous lesgo.conf
    def __init_folder(self, **kwargs):
        defaultKwargs = {
            'adjoint_flag'  : False,
            'dir'           : './'
        }
        kwargs = { **defaultKwargs, **kwargs }
        
        
        return

    def read_sensor_adjoint(self, t_ind):
        self.data['adjoint'] = np.empty((self.nx, self.nsensors, self.dims[0], self.dims[1], self.dims[2]))
        for nnx in range(self.nx):
            for ns in range(self.nsensors):
                self.data['adjoint'][nnx, ns, :, :, :] = read_array_from_file(self.adjoint_f % (ns+1, t_ind), self.dims)

        return
    
    def read_debug(self, t_ind):
        self.test_dir = self.output_dir + '/test'
        self.debug_advection_scalar(t_ind)
        self.debug_diffusion_scalar(t_ind)
        
        return

    def debug_advection_scalar(self, t_ind):
        adv_f = self.test_dir + r'/advection.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        
        dTdx_f = self.test_dir + r'/dTdx.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        dTdy_f = self.test_dir + r'/dTdy.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        dTdz_f = self.test_dir + r'/dTdz.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        
        self.data['adv'] = self._read_multi_scalar(adv_f, t_ind, self.dims)

        self.data['dTdx'] = self._read_multi_scalar(dTdx_f, t_ind, self.ihalf_dims)
        self.data['dTdy'] = self._read_multi_scalar(dTdy_f, t_ind, self.jhalf_dims)
        self.data['dTdz'] = self._read_multi_scalar(dTdz_f, t_ind, self.dims)

        return

    def debug_diffusion_scalar(self, t_ind):
        diff_f = self.test_dir + r'/diffusion.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        
        u_ihalf_f = self.test_dir + r'/u_ihalf.' + self.fmt_kwargs['fmt_tstep']
        v_jhalf_f = self.test_dir + r'/v_jhalf.' + self.fmt_kwargs['fmt_tstep']
        w_kw_f = self.test_dir + r'/w_kw.' + self.fmt_kwargs['fmt_tstep']
        
        theta_ihalf_f = self.test_dir + r'/theta_ihalf.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        theta_jhalf_f = self.test_dir + r'/theta_jhalf.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        theta_kw_f = self.test_dir + r'/theta_kw.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        
        self.data["diff"] = self._read_multi_scalar(diff_f, t_ind, self.dims)
        
        self.data["u_ihalf"] = read_array_from_file(u_ihalf_f % (t_ind), self.ihalf_dims)
        self.data["v_jhalf"] = read_array_from_file(v_jhalf_f % (t_ind), self.jhalf_dims)
        self.data["w_kw"] = read_array_from_file(w_kw_f % (t_ind), self.dims)
        
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
            'nk'    : 0,            # Default plot the obs of the first sensor
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
        axes[0].scatter(self.source_coords[self.nk, 0], self.source_coords[self.nk, 2], marker='x', color='red', alpha=0.7, )
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
            axes[i+1].plot(self.obs_t, self.data['obs'][i, kwargs['nk'], :], color='black', label='obs %.2i' % i, linewidth=0.5)
            axes[i+1].set_xlabel('t')
            axes[i+1].set_ylabel('theta')
            axes[i+1].set_ylim(kwargs['vmin']-0.1*abs(kwargs['vmax'] -  kwargs['vmin']), kwargs['vmax'])
            axes[i+1].set_xlim(self.obs_t.min(), self.obs_t.max())
            axes[i+1].legend(frameon=False)
            
        return fig, axes
    
    
    def sensor_measurements(self, sensor_locs, tt, **kwargs):
        import imageio
        defaultKwargs = {
            'gif_fname' : './result', 
            'gen_gif'   : True,
            'nk'        : 0
        }
        kwargs = { **defaultKwargs, **kwargs }
        
        self.nk = kwargs['nk'] # Assign value for the universal nk, to make sure that the obs of same sensor is used for reconstruction
        
        kwargs['gif_fname'] += ('_' + self.fmt_kwargs['fmt_ntheta'] + '_%s.gif' ) % ( kwargs['nk'], date_marker)
        
        self.sensor_init_(sensor_locs)
        self.obs_t = tt
        
        mesh = meshgrid(domain=self.domain, dims = self.dims)
        xi, yi, zi = mesh.cartesian2index(sensor_locs)
        
        if kwargs['gen_gif']:
            gif_dir = './tmp_gif'
            if not os.path.exists(gif_dir):
                os.makedirs(gif_dir)
            filenames = []
    
    
        self.data['obs'] = np.full(shape=(self.nsensors, self.ntheta, len(tt)), fill_value=np.nan)
        for t_ind, t in enumerate(tt):
            self.read_scalar(t) 
            for nT in range(self.ntheta):
                for ns in range(self.nsensors):
                    self.data['obs'][ns, nT, t_ind] = self.data['theta'][nT, xi[ns], yi[ns], zi[ns]] 
            
            if kwargs['gen_gif']:
                fig, ax = self.channel_obs(**kwargs)
                # create file name and append it to a list
                filename = gif_dir + f'/%.5i.png' % t
                filenames.append(filename)
                
                # save frame
                fig.savefig(filename, bbox_inches='tight')
                plt.close()
            
        if kwargs['gen_gif']:
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
        
    def sensor_field(self, **kwargs):
        defaultKwargs = {
            'axis_index'            : 0,
            'steps'                 : 16,
            'fieldname'             : "theta.IC",
            'readme_text'           : '',
            'input_dir'             : self.inputs_dir,
            'field_func'            : self.gaussian_field,
            'variance'              : 5e-2,
            'nx'                    : 1
            
        }
        kwargs = { **defaultKwargs, **kwargs }
        kwargs['readme_text'] += 'Domain =' + str(self.domain) + ', Dims = ' + str(self.dims) + ', \n Sensor source at %s' % str(self.s_locs)
        
        # self.xx = np.arange(0, len(self.coords[kwargs['axis_index']]), kwargs['steps'])
        # self.nx = self.xx.shape[0]
        self.nx = kwargs['nx']
        xsteps = int( self.dims[kwargs['axis_index']] / kwargs['nx'])

        self.xx = np.arange(0, len(self.coords[kwargs['axis_index']]), xsteps)
        
        # Original Sensor IC
        ic_data = [kwargs['field_func'](fieldname = kwargs['fieldname'], nk=i+1, source_point = [self.s_locs[0][i], self.s_locs[1][i], self.s_locs[2][i]], gen_fig=False, variance=kwargs['variance']) for i in range(self.nsensors)]

        rr = []
        for xind, xstep in enumerate(self.xx):
            for ind, data in enumerate(ic_data):
                shifted_data = np.roll(data, xstep, kwargs['axis_index'])
        
                fname = (kwargs['input_dir'] + '/%s.' + self.fmt_kwargs['fmt_ntheta']) %(kwargs['fieldname'], xind*len(ic_data) + ind + 1)
                
                rr.append(write_array_to_file(fname, shifted_data, **kwargs))
    
        return rr
    
    
    def adjoint_intersect(self, adjoint_tend = 1000, **kwargs):
        from itertools import combinations
        import imageio
        defaultKwargs = {
            'gif_fname' : './result',
        }
        kwargs = { **defaultKwargs, **kwargs }
        kwargs['gif_fname'] += ('_' + self.fmt_kwargs['fmt_ntheta'] + '_%s.gif' ) % ( self.nk, date_marker)
        
        
        
        adjoint_tind = int(np.argwhere(self.obs_t==adjoint_tend))
        recon_p = np.zeros(self.dims)
        
        gif_dir = './tmp_gif'
        if not os.path.exists(gif_dir):
            os.makedirs(gif_dir)
        filenames = []

        tind_max = self.obs_t.shape[0]
        for t_ind, t in enumerate(self.obs_t[:adjoint_tind]):
            self.read_sensor_adjoint(t)
            # for nnx in range(self.nx):
            self.nnx = 0
            for C_ind in combinations(range(self.nsensors), 2):
                    adjoint_ratio = self.data['adjoint'][self.nnx, C_ind[0]] / self.data['adjoint'][self.nnx, C_ind[1]]
                    obs_ratio = self.data['obs'][C_ind[0], self.nk, tind_max - t_ind - 1] /self.data['obs'][C_ind[1], self.nk, tind_max - t_ind - 1]
                    min_ratio = np.min( abs(adjoint_ratio - obs_ratio ))
                    recon_p += np.nan_to_num(min_ratio / abs(adjoint_ratio - obs_ratio ))
                    
                    fig, ax = self.channel_intersect(recon_p, t_ind = t_ind, **kwargs)
                    
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
            
        return 
    
    # [Working On]
    def channel_intersect(self, recon_p, **kwargs):
        from matplotlib import cm, colors, ticker
        import matplotlib.ticker as tick
        from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
        defaultKwargs = {
            'alpha' : 0.5,
            'vmin' : self.data['adjoint'].min(),
            'vmax' : self.data['adjoint'].max(),
            'cmap' : [cm.Reds, cm.Greens, cm.Blues],
            'levels': 101,
            'figsize' : (8, 8),
            'dpi'   : 150,
            'nk'    : 1,
            'y_ind' : 64,
            'clip'  : False,
            'norm'  : colors.Normalize,
            'antialiased' : False
        }
        kwargs = { **defaultKwargs, **kwargs }

        levels = np.linspace(kwargs['vmin'], kwargs['vmax'], kwargs['levels'])

        fig, axes = plt.subplots(3, 1, figsize=kwargs['figsize'], dpi = kwargs['dpi'])
        axes = axes.flatten()


        images = []
        [images.append(axes[0].contourf(self.coords[0], self.coords[2], self.data['adjoint'][self.nnx, i, :, kwargs['y_ind'], :].T, levels=kwargs['levels'],
                    norm=kwargs['norm'](vmin=kwargs['vmin'], vmax=kwargs['vmax'], clip=kwargs['clip']), cmap=kwargs['cmap'][i], alpha = kwargs['alpha'], antialiased=kwargs['antialiased']
                    )) for i in range(self.nsensors)]

        axes[0].scatter(self.source_coords[self.nk, 0], self.source_coords[self.nk, 2], marker='x', color='red', alpha=0.7, )
        axes[0].scatter(self.s_locs[0], self.s_locs[-1], marker='o', color='black', s = 1)
        axes[0].set_xlabel('x')
        axes[0].set_ylabel('z')
        axes[0].set_aspect('equal', 'box')
        axes[0].set_ylim(bottom=0)

        caxes = []
        for i in range(self.nsensors):
            ax1_divider = make_axes_locatable(axes[0])
            # Add an Axes to the right of the main Axes.
            caxes.append(ax1_divider.append_axes("right", size="2%", pad="1%"))
            cbar = fig.colorbar(
                cm.ScalarMappable(norm=kwargs['norm'](vmin=kwargs['vmin'], vmax=kwargs['vmax'], clip=kwargs['clip']), cmap=kwargs['cmap'][i]),
                extend='both', cax = caxes[i],
                                )

        # cbar.set_ticks(levels[::10])
        # cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))

        for i in range(self.nsensors):
            axes[1].plot(self.obs_t, self.data['obs'][i, self.nk, :], alpha=0.5, label='obs %.2i' % i, linewidth=0.5)
            # axes[1].set_yscale('log')
        axes[1].set_xlabel('t')
        axes[1].set_ylabel('theta')
        axes[1].set_ylim(kwargs['vmin']-0.1*abs(kwargs['vmax'] -  kwargs['vmin']), kwargs['vmax'])
        axes[1].set_xlim(self.obs_t.min(), self.obs_t.max())
        axes[1].legend(frameon=False)
        axes[1].vlines(self.obs_t[-kwargs['t_ind']-1], kwargs['vmin']-0.1*abs(kwargs['vmax'] -  kwargs['vmin']), kwargs['vmax'], color='red', linestyle =':')
            

        recon_x, recon_y, recon_z = np.unravel_index(recon_p.argmax(), self.dims)
        ca2 = axes[2].contourf(self.coords[0], self.coords[2], recon_p[:, kwargs['y_ind'], :].T, cmap = cm.Greys, alpha = 0.5, levels=100)
        axes[2].scatter(self.s_locs[0], self.s_locs[-1], marker='o', color='black', s = 1) # Sensor location
        axes[2].scatter(self.source_coords[self.nk, 0], self.source_coords[self.nk, 2], marker='x', color='red', alpha=0.7, )
        axes[2].set_xlabel('x')
        axes[2].set_ylabel('z')
        axes[2].set_aspect('equal', 'box')
        axes[2].set_ylim(bottom=0)
        axes[2].set_title('x=%.3f, y=%.3f, z=%.3f' %(self.coords[0][recon_x], self.coords[1][recon_y], self.coords[2][recon_z]))
        ax2_divider = make_axes_locatable(axes[2])
        # Add an Axes to the right of the main Axes.
        cax2 = ax2_divider.append_axes("right", size="2%", pad="1%")
        cbar = fig.colorbar(ca2, cax=cax2)
        return fig, axes
    
    
    
    # Pre-processing
    
    def install_lesgo(self, **repo_args):
        defaultKwargs = {
            'version'       : '6c4a959a5796127db5670f3357bc7608c67b2ba6',
            'repo'          : 'git@github.com:Advanced-Data-Assimilation/lesgo_eri.git',
            'dir'           : self.forward_dir + '/lesgo',
            'ssh_key'       : '/home/zyou6474/.ssh/id_rsa'
        }
        repo_args = { **defaultKwargs, **repo_args }
        self.repo_args = repo_args
        self.lesgo_dir = self.forward_dir + '/lesgo'
        
        
        import subprocess


        cmd = "unset SSH_ASKPASS; git clone %s %s; cd %s; git checkout %s" % (repo_args['repo'], repo_args['dir'], repo_args['dir'], repo_args['version'])

        returned_output = subprocess.check_output(cmd, shell=True)  # returns the exit code in unix
        # using decode() function to convert byte string to string
        print(returned_output.decode("utf-8"))
        
        cmd = "cd %s/src; ./build-lesgo" % (repo_args['dir'])

        build_output = subprocess.check_output(cmd, shell=True)  # returns the exit code in unix
        # using decode() function to convert byte string to string
        print(build_output.decode("utf-8"))
        
        cmd = "cp %s/build/lesgo-mpi-scalars-DA %s/" % (repo_args['dir'], self.forward_dir)

        returned_value = subprocess.run(cmd, shell=True)  # returns the exit code in unix
        # using decode() function to convert byte string to string
        print(returned_value)
        
        self.configs_dir = self.forward_dir + '/output/configs'
        if not os.path.exists(self.configs_dir):
            os.makedirs(self.configs_dir)
            print('create root folder in %s' % self.configs_dir)
            
            
        self._read_lesgoconf('%s/src/lesgo.conf' % repo_args['dir'])
        self._write_conf('%s/config_%s.ini' %(self.configs_dir, date_marker))        
        return
    
    def _read_lesgoconf(self, lesgoconf_f):
        lesgoconf = read_lesgoconf(lesgoconf_f)
        config = ConfigParser()
        for block in lesgoconf:
            print(block)
            config.add_section(block)
            for var in lesgoconf[block].keys():
                # print(var)
                config.set(block, var, lesgoconf[block][var])
                
        self.config = config
        print('Configs is loaded from %s' % lesgoconf_f)
        
        return
    
    def _write_conf(self, o_fname):
        configs_dir = os.path.dirname(o_fname)
        if not os.path.exists(configs_dir):
            os.makedirs(configs_dir)
            print('create root folder in %s' % configs_dir)
        with open(o_fname, 'w') as f:
            self.config.write(f)
        print('Configs is written in %s' % o_fname)
        
        return

    def _read_conf(self, conf):
        self.config = ConfigParser()
        self.config.read(conf)
        return

        
        
        
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
        self.source_coords = np.vstack((self.source_coords, np.array([x0, y0, z0])))
        
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
        
        return write_array_to_file(fname, data, **kwargs)
    
    def source_sameasic(self,):
        import shutil            
        ic_files = [filename for filename in os.listdir(self.inputs_dir) if filename.startswith("theta")]
        new_fname = ['/source.'+fname.split('.')[-1] for fname in ic_files]
        [shutil.copyfile(self.inputs_dir+ '/' + ic_files[ind], self.inputs_dir+name) for ind, name in enumerate(new_fname)]
        return new_fname
        
        

        
if __name__ == "__main__":
    lesgoconf_f = r'/home/zyou6474/Projects/lesgo_eri/src/lesgo.conf'
    

                
                    
                
