import numpy as np

str_func = lambda Lz, z, str_factor: Lz*(1+(np.tanh(str_factor*(z/Lz-1))/np.tanh(str_factor)))

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


class meshgrid():
    def __init__(self, domain, shape):
        
        self.domain = domain
        self.shape = shape
        
        
        self.coords = coords_xyz(domain, shape)
        self.dx = self.coords[0][1] - self.coords[0][0]

        
    def cartesian2index(self, points):
        px, py, pz = points
        
        x_ind = []
        y_ind = []
        z_ind = []
        
        for ind, _ in enumerate(px):
            dis_x = abs(self.coords[0] - px[ind])
            dis_y = abs(self.coords[1] - py[ind])
            dis_z = abs(self.coords[2] - pz[ind])
            
            if min(dis_x) + min(dis_y) + min(dis_z) < 3*self.dx:
                x_ind.append(np.argmin(dis_x))
                y_ind.append(np.argmin(dis_y))
                z_ind.append(np.argmin(dis_z))
            else:
                print('The %ith sensor locates outside mesh grid domain' % ind)
        
        self.points_ind = [x_ind, y_ind, z_ind]
        return x_ind, y_ind, z_ind
    
    
if __name__ == "__main__":
    mesh = meshgrid()
    sx = (4/3*np.pi, 4/3*np.pi, 4/3*np.pi)
    sy = (1/2*np.pi, 1/2*np.pi, 1/2*np.pi)
    sz = (1/4, 2/4, 3/4)
    points = (sx, sy, sz)
    
    xi, yi, zi = mesh._cartesian2index(points)