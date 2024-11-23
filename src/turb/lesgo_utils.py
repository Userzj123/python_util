import numpy as np
import matplotlib.pyplot as plt
import pyutils.plot_utils as pltutils
from pyutils.cartesian import coords_xyz, meshgrid, jaco
import os
from configparser import ConfigParser

from datetime import datetime
date_marker = datetime.today().strftime('%Y_%m_%d')

def read_periodic(f_half):
    half_dims = f_half.shape
    dims = (half_dims[0], half_dims[1], half_dims[2]*2)

    f = np.zeros(shape=dims)
    for i in range(half_dims[2]):
        f[:, :, i] = f_half[:, :, i]
        f[:, :, -i-1] = f_half[:, :, i]

    return f

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
        array = np.fromfile(f, dtype=dtype)

    array = array.reshape(shape, order='F')


    return array

def add_dirac_source(shape, ind_xset, ind_yset, ind_zset, volume, ratio=10):
    field = np.zeros(shape)
    nfield = len(ind_xset)
    for ind_field in range(nfield):
        # ldata.data['thetas_IC'][ind_source] = zz * 0
        sum_volume = 0
        for i in range(-2, 3):
            for j in range(-2, 3):
                for k in range(-2, 3):
                    location_ind = sorted([abs(i), abs(j), abs(k)])
                    if location_ind == [0, 0, 0]:
                        weight = 1
                    else:
                        weight = 0


                    field[ind_field, ind_xset[ind_field]+i, ind_yset[ind_field]+j,  ind_zset[ind_field]+k] = weight/volume[ind_xset[ind_field]+i, ind_yset[ind_field]+j,  ind_zset[ind_field]+k]
                    sum_volume += weight
        field /= sum_volume
    return field

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

def write_lesgoconf(root_dir, config: ConfigParser):


    for module in ['SCALAR_ADJOINT', 'SCALAR'] :
        if config.has_section(module):
            if module == 'SCALAR':
                fname  = root_dir + "/lesgo-scalar.conf"
            elif module == "SCALAR_ADJOINT":
                fname  = root_dir + "/lesgo-scalar-adjoint.conf"

            with open(fname, 'w') as f:
                for sec in [module, module+"_LOG"]:
                    f.write(sec + '{ \n')
                    for var in config[sec]:
                        f.write(var + ' = ' + config[sec][var] + '\n')
                        
                    f.write('} \n\n')

                    config.remove_section(sec)
            print('Write configs in %s' % fname)            

    fname  = root_dir + "/lesgo.conf"
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
    
    def __init__(self, domain, dims, root_dir, ntheta=1) -> None:
        self.root_dir = root_dir
        if not os.path.exists(root_dir):
            os.makedirs(root_dir)
            print('create root folder in %s' % root_dir)
        self.inputs_dir = root_dir + '/inputs'
        if not os.path.exists(self.inputs_dir):
            os.makedirs(self.inputs_dir)
            print('create inputs folder in %s' % self.inputs_dir)
            


        self.output_dir = root_dir + '/output'
        
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
        
        iijjkk_dims = dims*1
        iijjkk_dims = [i+1 for i in iijjkk_dims]

        coords = coords_xyz(domain, iijjkk_dims, center=False, stretch=True)
        
        self.volume = np.zeros(shape=dims)
        z = np.linspace(0, domain[2], dims[2]+1)
        z = (z[1:] + z[:-1])/2
        dz_uv = jaco(domain[2], z) * (z[1] - z[0])
        for i in range(dims[0]):
            for j in range(dims[1]):
                for k in range(dims[2]):
                    self.volume[i, j, k] = (coords[0][i+1] - coords[0][i]) * (coords[1][j+1] - coords[1][j]) * dz_uv[k]
        

        # for i in range(dims[0]):
        #     for j in range(dims[1]):
        #         for k in range(dims[2]):
        #             self.volume[i, j, k] = (coords[0][i+1] - coords[0][i]) * (coords[1][j+1] - coords[1][j]) * (coords[2][k+1] - coords[2][k])
        
        self.ihalf_coords = coords_xyz(domain, self.ihalf_dims, center=False, stretch=True)
        self.jhalf_coords = coords_xyz(domain, self.jhalf_dims, center=False, stretch=True)
        self.kw_coords = coords_xyz(domain, self.kw_dims, center=False, stretch=True)
        
        self._fnames()
        return
    # Post-processing
    
    
    def read_data(self, t_ind):
        if not self.adjoint:
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
    
    def read_scalar(self, t_ind):
        thetas = read_array_from_file(self.theta_f % t_ind, self.source_dims)
        
        self.data['theta'] = thetas
        return thetas

    def read_scalar_adjoint(self, t_ind):
        theta_Ds = read_array_from_file(self.adjoint_f % t_ind, self.sensor_dims)
        
        self.data['theta_D'] = theta_Ds
        return
    
    def read_inputs(self, dir, in_dims = None):
        """
        Read the data of inputs field.
        """
        self.read_velocity_inputs(dir, in_dims = in_dims, in_domain=self.domain)
        
        self.read_scalar_inputs(dir, )
        

        return

    def read_scalar_inputs(self, dir):
        source_fs = [os.path.join(dir, filename) for filename in os.listdir(dir) if filename.startswith("source.")]
        theta_IC_fs = [os.path.join(dir, filename) for filename in os.listdir(dir) if filename.startswith("theta.")]

        self.data['source'] = np.empty((0, self.dims[0], self.dims[1], self.dims[2]))
        for find, fn in enumerate(source_fs):
            self.data['source'] = np.concatenate((self.data['source'], read_array_from_file(fn, self.dims)[np.newaxis, :, :, :]))
        
        
        self.data['theta_IC'] = np.empty((0, self.dims[0], self.dims[1], self.dims[2]))
        for find, fn in enumerate(theta_IC_fs):
            self.data['theta_IC'] = np.concatenate((self.data['theta_IC'], read_array_from_file(fn, self.dims)[np.newaxis, :, :, :]))
        return
    



    def read_velocity_inputs(self, dir, in_dims = None, in_domain=None, periodic=False):
        from scipy.interpolate import interpn
        """_summary_

        Args:
            dir (_type_): folder that store the initial condition.
            in_dims (_type_, optional): Dimension of read in initial condition. Defaults to None.
        """

        
        u_icf = dir + '/u_velocity.IC'
        v_icf = dir + '/v_velocity.IC'
        w_icf = dir + '/w_velocity.IC'

        if periodic:
            half_dims = self.dims * 1
            half_dims[2] = int(self.dims[2]/2)


            u_ic_half = read_array_from_file(u_icf, half_dims)
            v_ic_half = read_array_from_file(v_icf, half_dims)
            w_ic_half = read_array_from_file(w_icf, half_dims)

            self.data['u_ic'] = read_periodic(u_ic_half)
            self.data['v_ic'] = read_periodic(v_ic_half)
            self.data['w_ic'] = read_periodic(w_ic_half)
        else:

            if in_dims == None:
                self.data['u_ic'] = read_array_from_file(u_icf, self.dims)
                self.data['v_ic'] = read_array_from_file(v_icf, self.dims)
                self.data['w_ic'] = read_array_from_file(w_icf, self.dims)
            else:
                u_ic = read_array_from_file(u_icf, in_dims)
                v_ic = read_array_from_file(v_icf, in_dims)
                w_ic = read_array_from_file(w_icf, in_dims)
                
                in_coords = coords_xyz(in_domain, in_dims, center=False, stretch=True)
                points = np.array([np.meshgrid(*self.coords, indexing="ij")[ind].reshape(-1) for ind in range(3)])
                
                self.data['u_ic'] = interpn(in_coords, u_ic, points.T).reshape(*self.dims)
                self.data['v_ic'] = interpn(in_coords, v_ic, points.T).reshape(*self.dims)
                self.data['w_ic'] = interpn(in_coords, w_ic, points.T).reshape(*self.dims)
            
        return 1
            
    
    def write_inputs(self):
        # Velocity Initial Conditions
        write_array_to_file(self.inputs_dir + '/u_velocity.IC', self.data['u_ic'])
        write_array_to_file(self.inputs_dir + '/v_velocity.IC', self.data['v_ic'])
        write_array_to_file(self.inputs_dir + '/w_velocity.IC', self.data['w_ic'])
        
        # Scalar Initial Conditions
        self.nsource = self.data['thetas_IC'].shape[0]
        self.source_dims = self.dims*1
        self.source_dims.insert(0, self.nsource)
        self.source_dims = tuple(self.source_dims)
        self.config['SCALAR']['ntheta'] = (self.fmt_kwargs['fmt_ntheta']) % self.nsource
        thetaic_f = self.inputs_dir + '/thetas.IC'
        write_array_to_file(thetaic_f, self.data['thetas_IC'])
            
        # Source
        if self.config['SCALAR']['opt_source'] == '.TRUE.':
            if self.data['sources'].shape[0] != self.nsource:
                raise Exception("Sources of all scalar need to be defined in self.data['source'] with shape of (nk, nx, ny, nz)")
            else:   
                source_f = self.inputs_dir + '/sources'
                write_array_to_file(source_f, self.data['sources'])

        # Sensor
        if self.config.has_section('SCALAR_ADJOINT'):
            self.nsensor = self.data['theta_Ds_IC'].shape[0]
            self.sensor_dims = self.dims*1
            self.sensor_dims.insert(0, self.nsensor)
            self.sensor_dims = tuple(self.sensor_dims)
            self.config['SCALAR_ADJOINT']['nadjoint'] = (self.fmt_kwargs['fmt_ntheta']) % self.nsensor

            sensor_f = self.inputs_dir + '/theta_Ds.IC'
            write_array_to_file(sensor_f, self.data['theta_Ds_IC'])

            if self.config['SCALAR_ADJOINT']['opt_source_d'] == '.TRUE.':
                sensor_source_f = self.inputs_dir + '/source_Ds'
                write_array_to_file(sensor_source_f, self.data['source_Ds'])
            
        # LESGO.CONF
        ## Update shape of domain 
        self.config['DOMAIN']['nx'] = f'{self.dims[0]}'
        self.config['DOMAIN']['ny'] = f'{self.dims[1]}'
        self.config['DOMAIN']['nz'] = f'{self.dims[2]}'
        self.config['DOMAIN']['lx'] = f'{self.domain[0]}'
        self.config['DOMAIN']['ly'] = f'{self.domain[1]}'
        self.config['DOMAIN']['lz'] = f'{self.domain[2]}'
        ## Write lesgo.conf
        write_lesgoconf(self.inputs_dir, self.config)
        
        return 

    def _fnames(self, **fmt_kwargs):
        defaultKwargs = {
            'fmt_tstep'   : r'%.8i',
            'fmt_ntheta'   : r'%.3i'
        }
        self.fmt_kwargs = { **defaultKwargs, **fmt_kwargs }
        

        self.u_f = self.output_dir + r'/ns/u_velocity.' + self.fmt_kwargs['fmt_tstep']
        self.v_f = self.output_dir + r'/ns/v_velocity.' + self.fmt_kwargs['fmt_tstep']
        self.w_f = self.output_dir + r'/ns/w_velocity.' + self.fmt_kwargs['fmt_tstep']

        
        self.theta_f = self.output_dir + r'/scalar/thetas.'+ self.fmt_kwargs['fmt_tstep']
        self.adjoint_f = self.output_dir + r'/scalar_adjoint/theta_Ds.' + self.fmt_kwargs['fmt_tstep']
            
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

        
        

        
if __name__ == "__main__":
    root_dir = '/home/zyou6474/tasks/source_inversion/forward'

    dims = [128, 128, 64]
    domain = [2*np.pi, np.pi, 1]
    ntheta = 1

    ldata = lesgo_data(domain, dims, root_dir, ntheta=ntheta)
    ldata._fnames(fmt_ntheta='%.3i')

    t_total = 500
    tstep = 1
    tt = np.arange(0, t_total+1, tstep)

    sensor_locs = np.array([
        [1.6 * np.pi, 1.6 * np.pi, 1.6 * np.pi,],
        [0.5 * np.pi, 0.5 * np.pi, 0.5 * np.pi,],
        [0.3, 0.5, 0.7]
    ])

    from pyutils.cartesian import coords_xyz, meshgrid
    mesh = meshgrid(domain=ldata.domain, shape = ldata.dims)
    xi, yi, zi = mesh.cartesian2index(sensor_locs)

    ldata.sensor_measurements(sensor_locs, tt)
    

                
                    
                
