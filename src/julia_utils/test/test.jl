include("../lst.jl")
using .lst_scalar

using LinearAlgebra
using Plots
using FFTW
using MAT
plotlyjs()


lst = generate_LST(384, 2)

fname = "/Users/user/Documents/Projects/python_util/data/scalar_LST_matlab_linux.mat"
mat_vars = matread(fname)

scatter(mat_vars["omega"])
scatter!(1im .* lst.F.values)