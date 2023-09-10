module mesh

export coordinates, init_mesh, meshgrid

struct coordinates
    x :: Vector
    y :: Vector
    z :: Vector
end

struct surface
    Sx :: Array
    Sy :: Array
    Sz :: Array
end

struct meshgrid
    ijk :: coordinates
    ijkk :: coordinates

    S :: surface
    V :: Array
end


function generate_coords(domain, dims; center=false, stretch=false, str_factor=1.5)

    if center
        # Generate 1D arrays of x, y, and z coordinates
        x_coords = LinRange(0, domain[1], dims[1]+1)
        y_coords = LinRange(0, domain[2], dims[2]+1)
        z_coords = LinRange(0, domain[3], dims[3]+1)

        x_coords = 0.5 * x_coords[1:end-1] + 0.5 * x_coords[2:end]
        y_coords = 0.5 * y_coords[1:end-1] + 0.5 * y_coords[2:end]
        z_coords = 0.5 * z_coords[1:end-1] + 0.5 * z_coords[2:end]
        
    else
        # Generate 1D arrays of x, y, and z coordinates
        x_coords = LinRange(0, domain[1], dims[1])
        y_coords = LinRange(0, domain[2], dims[2])
        z_coords = LinRange(0, domain[3], dims[3])
    end

    # k(uv) grid
    if stretch
        z_stretch = domain[3] .* (1 .+(tanh.(str_factor .*(z_coords ./domain[3] .-1)) ./tanh(str_factor)))

        coords = coordinates(x_coords, y_coords, z_stretch)
        return coords
    else
        coords = coordinates(x_coords, y_coords, z_coords)
        return coords
    end

end


function generate_meshgrid(domain, dims)::meshgrid
end

end