module topology

export coordinates, generate_coords, generate_meshgrid

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
    xx :: Array{Float64, 3}
    yy :: Array{Float64, 3}
    zz :: Array{Float64, 3}
end

struct mdata
    ijk :: meshgrid
    ijkk :: meshgrid

    S :: surface
    V :: Array
end


function generate_coords(domain, dims; center=false, stretch=false, str_factor=1.5)::coordinates

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


function generate_meshgrid(coord::coordinates)::meshgrid
    x = reshape(coord.x, (length(coord.x), 1, 1))
    y = reshape(coord.y, (1, length(coord.y), 1))
    z = reshape(coord.z, (1, 1, length(coord.z)))

    xx = repeat(x, 1, length(y), length(z))
    yy = repeat(y, length(x), 1, length(z))
    zz = repeat(z, length(x), length(y), 1)
    return meshgrid(xx, yy, zz)
end

end