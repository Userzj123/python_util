include("./mesh.jl")
using .mesh


domain = [2*pi, pi, 1]
dims = [16, 16, 128]

coord= init_mesh(domain, dims)
