module lesgo


using LinearAlgebra

struct point_index
    x::Int16
    y::Int16
    z::Int16
end

struct points_index
    points::Array{point_index}
end


# struct coord

# end

function coords_xyz(domain, dims, stretch=false, str_factor=1.5, center=false)
    # Define the range and number of points in each direction
    x_range = (0, domain[1])
    y_range = (0, domain[2])
    z_range = (0, domain[3])
    n_x, n_y, n_z = dims

    if center
        # Generate 1D arrays of x, y, and z coordinates
        x_coords = LinRange(x_range[1], x_range[2], Int(n_x)+1)
        y_coords = LinRange(y_range[1], y_range[2], Int(n_y)+1)
        z_coords = LinRange(z_range[1], z_range[2], Int(n_z)+1)

        x_coords = 0.5 * x_coords[1:end-1] + 0.5 * x_coords[2:end]
        y_coords = 0.5 * y_coords[1:end-1] + 0.5 * y_coords[2:end]
        z_coords = 0.5 * z_coords[1:end-1] + 0.5 * z_coords[2:end]
        
    else
        # Generate 1D arrays of x, y, and z coordinates
        x_coords = LinRange(x_range[1], x_range[2], Int(n_x))
        y_coords = LinRange(y_range[1], y_range[2], Int(n_y))
        z_coords = LinRange(z_range[1], z_range[2], Int(n_z))
    end

    # k(uv) grid
    if stretch
        z_stretch = z_range[2]*(1 .+(tanh.(str_factor*(z_coords/z_range[2] .-1))/tanh.(str_factor)))
        return x_coords, y_coords, z_stretch
    end 

    return x_coords, y_coords, z_coords
end 

end