include("./mesh.jl")
using .mesh
using TOML

domain = [2*pi, pi, 1]
dims = [16, 16, 128]

coord= init_mesh(domain, dims)


config = TOML.parsefile("/home/zejiany/Documents/Projects/python_util/src/julia_cfd/data/cfd.conf")