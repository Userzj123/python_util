module topology

export coordinates, generate_coords, generate_meshgrid

struct coordinates
    x :: Vector
    y :: Vector
    z :: Vector
    jaco_z :: Vector
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
        tmp_x = LinRange(0, domain[1], dims[1]+1)
        tmp_y = LinRange(0, domain[2], dims[2]+1)
        tmp_z = LinRange(0, domain[3], dims[3]+1)

        x_coords = 0.5 * tmp_x[1:end-1] + 0.5 * tmp_x[2:end]
        y_coords = 0.5 * tmp_y[1:end-1] + 0.5 * tmp_y[2:end]
        z_coords = 0.5 * tmp_z[1:end-1] + 0.5 * tmp_z[2:end]
        
    else
        # Generate 1D arrays of x, y, and z coordinates
        x_coords = LinRange(0, domain[1], dims[1]+1)
        y_coords = LinRange(0, domain[2], dims[2]+1)
        z_coords = LinRange(0, domain[3], dims[3]+1)
    end

    # k(uv) grid
    if stretch
        z_stretch = domain[3] .* (1 .+(tanh.(str_factor .*(z_coords ./domain[3] .-1)) ./tanh(str_factor)))
        # FIELD2(:) = L_z*(1.0_rprec+(tanh(str_factor*(z_uv(:)/L_z-1.0_rprec))    &
        # /tanh(str_factor)))
        jaco_z = domain[3] * (str_factor /domain[3]) .* (1 .- (tanh.(str_factor.*(z_coords./domain[3] .- 1))).^2) ./ tanh.(str_factor)
        # FIELD4(:) = L_z*(str_factor/L_z)*                                       &
        # (1-(tanh(str_factor*(z_uv(:)/L_z-1.0_rprec)))**2)/tanh(str_factor)

        coords = coordinates(x_coords, y_coords, z_stretch, jaco_z)
        return coords
    else
        jaco_z = z_coords .* 0 .+ 1
        coords = coordinates(x_coords, y_coords, z_coords, jaco_z)
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

function get_volume(coord::coordinates)::Array{Float64, 3}
    Nx = length(coord.x)
    Ny = length(coord.y)
    Nz = length(coord.z)


    volume = Array{Float64, 3}(undef, Nx-1, Ny-1, Nz-1)
    for i in range(1, Nx-1)
        for j in range(1, Ny-1)
            for k in range(1, Nz-1)
                volume[i, j, k] = (coord.x[i+1] - coord.x[i]) * (coord.y[j+1] - coord.y[j]) * (coord.z[k+1] - coord.z[k])
            end
        end
    end

    return volume
end

end