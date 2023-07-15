import numpy as np
import matplotlib.pyplot as plt
import pyutils.plot_utils as pltutils
from pyutils.cartesian import coords_xyz, meshgrid
import os

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

class lesgo_data():
    """
    Read lesgo's simulation result files.
    
    """
    
    def __init__(self, domain, dims, forward_dir, ntheta=1) -> None:
        if not os.path.exists(forward_dir):
            os.makedirs(forward_dir)
            print('create root folder in %s' % forward_dir)
        self.forward_inputs_dir = forward_dir + '/inputs'
        self.inputs_dir = self.forward_inputs_dir
        if not os.path.exists(self.inputs_dir):
            os.makedirs(self.inputs_dir)
            print('create inputs folder in %s' % self.inputs_dir)
            
        
        self.forward_outputs_dir = forward_dir + '/outputs'
        self.outputs_dir = self.forward_outputs_dir
        
        self.coords = coords_xyz(domain, dims, center=True, stretch=True)
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
    
    def lesgo_conf(self, **lesgo_args):
        default_lesgoargs = {
            'u_star'    : 1,
            'z_i'       : 1,
            'nu_molec'  : 0.00556,
            'default_conf' : '' #???
        }
        lesgo_args = { **default_lesgoargs, **lesgo_args }
        
        
        with open(lesgo_args['default_conf'], 'r') as f:
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
                        values = value.replace('\n', '').strip().split()
                        if len(values) == 1:
                            if values[0] == '.true.':
                                block_tmp[varname] = True
                            elif values[0] == '.false.':
                                block_tmp[varname] = False
                            else:
                                try: block_tmp[varname] = float(values[0])
                                except:  block_tmp[varname] = value
                        else:
                            try: block_tmp[varname] = [float(d) for d in values]
                            except:  block_tmp[varname] = value
                        
                    if '{' in line:
                        contents = line.split('{')
                        block_tmp['blockname'] = contents[0].strip()
                        block_entry = 1
                    elif '}' in line:
                        configs[block_tmp['blockname']] = block_tmp
                        block_entry = 2
                        block_tmp = {}
                    
                    
        return configs
    
    
    def read_data(self, t_ind):
        self.read_velocity(t_ind)
        self.read_scalar(t_ind)
        
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
    
    def read_inputs(self, ):
        """
        Read the data of inputs field.
        """
        if not self.adjoint :
            self.u_icf = self.inputs_dir + '/u_velocity.IC'
            self.v_icf = self.inputs_dir + '/v_velocity.IC'
            self.w_icf = self.inputs_dir + '/w_velocity.IC'
            
            self.data['u_ic'] = read_array_from_file(self.u_icf, self.dims)
            self.data['v_ic'] = read_array_from_file(self.v_icf, self.dims)
            self.data['w_ic'] = read_array_from_file(self.w_icf, self.dims)
        
        self.source_fs = [os.path.join(self.inputs_dir, filename) for filename in os.listdir(self.inputs_dir) if filename.startswith("source.")]
        self.theta_IC_fs = [os.path.join(self.inputs_dir, filename) for filename in os.listdir(self.inputs_dir) if filename.startswith("theta.")]

        self.data['source'] = np.empty((0, self.dims[0], self.dims[1], self.dims[2]))
        for find, fn in enumerate(self.source_fs):
            self.data['source'] = np.concatenate((self.data['source'], read_array_from_file(fn, self.dims)[np.newaxis, :, :, :]))
        
        
        self.data['theta.IC'] = np.empty((0, self.dims[0], self.dims[1], self.dims[2]))
        for find, fn in enumerate(self.theta_IC_fs):
            self.data['theta.IC'] = np.concatenate((self.data['theta.IC'], read_array_from_file(fn, self.dims)[np.newaxis, :, :, :]))
        
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
    
        if self.adjoint:
            self.adjoint_f = self.adjoint_dir + r'/outputs/theta.' + self.fmt_kwargs['fmt_ntheta'] + '.' + self.fmt_kwargs['fmt_tstep']
            
        return 
    
    
    def set_adjoint(self, adjoint = False, adjoint_dir = None, **kwargs):
        self.adjoint = adjoint
        if adjoint:
            self.adjoint_dir = adjoint_dir
            self.adjoint_inputs_dir = self.adjoint_dir + '/inputs'
            self.adjoint_outputs_dir = self.adjoint_dir + '/outputs'
            
            self.inputs_dir = self.adjoint_inputs_dir
            self.outputs_dir = self.adjoint_outputs_dir
            if not os.path.exists(adjoint_dir):
                os.makedirs(adjoint_dir)
                print('create adjoint root folder in %s' % adjoint_dir)
            if not os.path.exists(self.adjoint_inputs_dir):
                os.makedirs(self.adjoint_inputs_dir)
                print('create adjoint inputs folder in %s' % self.adjoint_inputs_dir)
        else:
            self.inputs_dir = self.forward_inputs_dir
            self.outputs_dir = self.forward_outputs_dir
                
        self._fnames(**kwargs)
        return

    def read_sensor_adjoint(self, t_ind):
        self.data['adjoint'] = np.empty((self.nx, self.nsensors, self.dims[0], self.dims[1], self.dims[2]))
        for nnx in range(self.nx):
            for ns in range(self.nsensors):
                self.data['adjoint'][nnx, ns, :, :, :] = read_array_from_file(self.adjoint_f % (ns+1, t_ind), self.dims)

        return
    
    def read_debug(self, t_ind):
        self.debug_advection_scalar(t_ind)
        self.debug_diffusion_scalar(t_ind)
        
        return

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
        
        u_ihalf_f = self.outputs_dir + r'/u_ihalf.' + self.fmt_kwargs['fmt_tstep']
        v_jhalf_f = self.outputs_dir + r'/v_jhalf.' + self.fmt_kwargs['fmt_tstep']
        w_kw_f = self.outputs_dir + r'/w_kw.' + self.fmt_kwargs['fmt_tstep']
        
        theta_ihalf_f = self.outputs_dir + r'/theta_ihalf.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        theta_jhalf_f = self.outputs_dir + r'/theta_jhalf.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        theta_kw_f = self.outputs_dir + r'/theta_kw.' + self.fmt_kwargs['fmt_ntheta'] +'.' + self.fmt_kwargs['fmt_tstep']
        
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
            'fieldname'             : "source",
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
    

                
                    
                
