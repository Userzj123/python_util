# include("../mesh.jl")
# using .topology
using LinearAlgebra
using FFTW

const domain = [2*pi, pi, 1]
const dims = [16, 16, 128]
const Pr = 0.71
const nu = 1/180/Pr
const dt = 1e-5


# coord= generate_coords(domain, dims)
# mesh = generate_meshgrid(coord)


const Nx = dims[1]
const Nz = dims[3]
const Ntheta = dims[3]

const dx = domain[1]/(dims[1]-1)
const dz = domain[3]/(dims[3]-1)

# const varstype = Array{Float64, 2}(undef, Nx, Nz)

# x = reshape(coord.x, (Nx, 1))
# z = reshape(coord.z, (1, Nz))

x = LinRange(0, domain[1], Nx)
z = LinRange(0, domain[3], Nz)



xx = repeat(x, 1, Nz)
zz = repeat(z', Nx, 1)

u = zz .* (2 .- zz)
v = zz .* 0
w = zz .* 0


theta_ic = zeros(Float64, (Nx, Nz, Ntheta))
[theta_ic[:, k, k] = sin.(x) for k in range(1, Nz)]


function solver(u, v, w, T)
    T[:, 1] .= 0
    T[:, end] .= 0

    kx = 1
    ky = 0

    # input = [fft(T[:, k], 1)[kx] for k in range(1, Nz)]
    input = fft(T, 2)[kx+1, :]


    RHS = zeros(Float64, (Nx, Nz))
    RHS[:, 2:end-1] .= T[:, 2:end-1] .- 
    dt * nu .*(2 .*T[:, 2:end-1] .- circshift(T[:, 2:end-1],[1,0]) .- circshift(T[:, 2:end-1],[-1,0]))./dx./dx .- 
    dt * nu .*(2 .*T[:, 2:end-1] .- T[:, 1:end-2] .- T[:, 3:end]) ./dz./dz .+ 
    dt * u[:, 2:end-1].*(circshift(T[:, 2:end-1],[1, 0]) .- circshift(T[:, 2:end-1],[-1, 0]))./dx./2

    output = fft(RHS, 2)[kx+1, :]
    return RHS, input, output
end



result = zeros(Float64, (Nx, Nz, Ntheta))
A = Array{ComplexF64, 2}(undef, Nz, Nz)
B = Array{ComplexF64, 2}(undef, Nz, Nz)

for k in range(1, Nz)
    result[:, :, k], A[:, k], B[:, k] = solver(u, v, w, theta_ic[:, :, k])
end

# [result[:, :, k] .= solver!(u, v, w, theta_ic[:, :, k]) for k in range(1, Nz)]

