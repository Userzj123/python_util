include("./lesgo_utils.jl")
import .lesgo
using MAT

domain = [2*pi, pi, 1]
dims = [16, 16, 128]

x, y, z = lesgo.coords_xyz(domain, dims, true, 1.5, true)

vars = matread("src/julia_utils/data/lesgo_coords.mat")
