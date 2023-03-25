import numpy as np 
import matplotlib.pyplot as plt

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

def gaussian_ic(domain, dims, mu_L=np.array([0.5, 0.5, 0.2]), sigma=np.eye(3)*1e-3):
    from scipy.stats import multivariate_normal
    X, Y, Z = XYZ(domain, dims)

    xyz = np.column_stack([X.flat, Y.flat, Z.flat])

    data = multivariate_normal.pdf(xyz, mean=mu_L*domain, cov=sigma)

    data = data.reshape(X.shape)
    return data


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

def velocity_parabolic(domain, dims):
    X, Y, Z = XYZ(domain, dims)

    u = Z**2
    v = X*0
    w = X*0

    return u, v, w


def dirac_source(domain, dims, source_z):
    Lx, Ly, Lz = domain
    nx, ny, nz = dims

    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    z = np.linspace(0, Lz, nz)

    Z, Y, X = np.meshgrid(z, y, x, indexing='ij')

    theta = X*0
    theta[source_z, 64, 64] += 1
    return theta

def multi_theta(domain, dims, ntheta=3):
    Lx, Ly, Lz = domain
    nx, ny, nz = dims

    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    z = np.linspace(0, Lz, nz)

    Z, Y, X = np.meshgrid(z, y, x,  indexing='ij')

    theta = np.zeros((ntheta, nz, ny, nx))
    theta[0, 20, 64, 64] += 1
    theta[1, 32, 64, 64] += 1
    theta[2, 50, 64, 64] += 1
    return theta

def theta_IC(domain, dims, path, ic_type="dirac_source", source_z=32):
    print(ic_type)

    if ic_type == "dirac_source":
        write_array_to_file(path % 'theta.IC', dirac_source(domain, dims, source_z))
        print('updated theta IC of forward simulation')
    elif ic_type == "gaussian_forward":
        write_array_to_file(path % 'theta.00000000', gaussian_ic(domain, dims, mu_L=np.array([0.5, 0.5, 0.2])))
        print('updated theta IC of forward simulation')
    elif ic_type == "gaussian_backward":
        write_array_to_file(path % 'theta.IC', gaussian_ic(domain, dims, mu_L=np.array([0.5, 0.5, 0.8])))
        print('updated theta IC of backward simulation')
    if ic_type == "multi_source":
        write_array_to_file(path % 'thetas.IC', multi_theta(domain, dims, 3))
        print('updated theta IC of forward simulation')
    if ic_type == "multi_source_seperate":
        write_array_to_file(path % 'theta.IC.01', dirac_source(domain, dims, 32))
        write_array_to_file(path % 'theta.IC.02', dirac_source(domain, dims, 20))
        write_array_to_file(path % 'theta.IC.03', dirac_source(domain, dims, 50))

        print('updated theta IC of forward simulation')

        

def velocity_IC(domain, dims, path, ic_type=None, base_path=None, base_dim=None):
    if ic_type == None: 
        print('Specify IC type')
        return 0
    
    if ic_type == "parabolic":
        u, v, w = velocity_parabolic(domain, dims)

        write_array_to_file(path % 'u_velocity.IC', u)
        write_array_to_file(path % 'v_velocity.IC', v)
        write_array_to_file(path % 'w_velocity.IC', w)
        print('updated velocity IC of forward simulation')

    if ic_type == "interpolation":
        interpolate_velocity_IC(domain, dims, path % 'u_velocity.00000000', base_dim, base_path % 'u_velocity.IC' )
        interpolate_velocity_IC(domain, dims, path % 'v_velocity.00000000', base_dim, base_path % 'v_velocity.IC' )
        interpolate_velocity_IC(domain, dims, path % 'w_velocity.00000000', base_dim, base_path % 'w_velocity.IC' )

        print('updated velocity IC of forward simulation')

def interpolate_velocity_IC(domain, out_dims, out_path, dims, path):
    from scipy.interpolate import RegularGridInterpolator as rgi

    if path == None or dims == None: print("Require base data")

    Lx, Ly, Lz = domain
    nx, ny, nz = dims

    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    z = np.linspace(0, Lz, nz)

    nX, nY, nZ = out_dims

    X, Y, Z = XYZ(domain, out_dims)

    with open(path, 'rb') as f:
        U_origin = np.fromfile(f, dtype=np.float64,)

    U_origin = np.reshape(U_origin, dims[::-1])


    my_interpolating_function = rgi((z, y, x), U_origin)
    U = my_interpolating_function(np.array([Z.reshape(-1), Y.reshape(-1), X.reshape(-1)]).T)

    # U = U.reshape(out_dims[::-1])
    # write_binary(out_path, U)
    # Previous version
    U.tofile(out_path)
    print("sucessfully upload interpolated velocity IC")
    return U

if __name__ == "__main__":

    dims = [128, 128, 64]
    domain = [2*np.pi, np.pi, 1]
    path = '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2_channel_flow_scalar/inputs/%s'

    out_dims = [256, 256, 128]
    out_path = '/Users/user/Documents/Projects/python_util/inputs/%s'

    
    # velocity_IC(domain, dims=out_dims, path=out_path, ic_type="interpolation", base_path=path, base_dim=dims)
    theta_IC(domain, out_dims, out_path, "gaussian_forward", source_z=64)
