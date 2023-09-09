module mesh

export coords, init_mesh

struct coordinates
    x :: Vector
    y :: Vector
    z :: Vector
end

global xx::Array, yy::Array, zz::Array
global coords::coordinates


function init_mesh(domain, dims)
    x = LinRange(0, domain[1], dims[1])
    y = LinRange(0, domain[2], dims[2])
    z = LinRange(0, domain[3], dims[3])


    coords = coordinates(x, y, z)


    return coords
end

end