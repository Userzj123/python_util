include("./topology.jl")
using .topology
using LinearAlgebra
using FFTW

const Pr = 0.71
const nu = 1/180/Pr
const dt = 1e-4


const Nx = 400
const Nz = 129

const kx = 1
const ky = 0


const domain = [2*pi, pi, 2]
const dims = [Nx, 16, Nz-2]
coords = topology.generate_coords(domain, dims; center=true, stretch=false)

x = coords.x
z = [-1 .* coords.z[1]; coords.z; 2*domain[3]-coords.z[end]]


const dx = x[2] - x[1]
const dz = z[2] - z[1]
const Lx = domain[1]

xx = x .+ z'.*0
zz = x .* 0 .+ z'

u = zz .* ( 2 .- zz)
v = zz .* 0
w = zz .* 0;


function implicit_solver(u, v, w, T)
    # Modified Wavenumber k = sin(kx * dx/2) * 2 /dx    
    mkx = zeros(Float64, (Nx));
    for ki in range(1, Nx)
        kkx = 2*pi/Lx*(ki-1);
        mkx[ki] = 2*abs(sin(0.5*kkx*dx))/dx;
    end


    RHS = zeros(Float64, (Nx, Nz));
    RHS .= T - dt .* u .* (circshift(T,[-1, 0]) .- circshift(T,[1, 0]))./dx./2;

    RHS_hat = zeros(ComplexF64, (Nx, Nz));
    for k in range(1, Nz)
        RHS_hat[:, k] = fft(RHS[:, k]);
    end

    D2 = diagm(0 => -2*ones((Nz)), 1 => ones((Nz-1)), -1 => ones((Nz-1)));
    D2 .= D2/dz/dz;


    T_hat = zeros(ComplexF64, (Nx, Nz));
    for ki in range(1, Nx)
        D = diagm(ones((Nz))) - dt*nu * (D2 - diagm(ones((Nz))) * mkx[ki]^2 );
        D[1, 1] = 1;
        D[1, 2] = 1;
        D[end, end-1] = 1;
        D[end, end] = 1;

        T_hat[ki, :] = D\RHS_hat[ki, :];
    end
    

    for k in range(1, Nz)
        T[:, k] = real(ifft(T_hat[:, k]));
    end

    return T
end


A = Array{ComplexF64, 2}(undef, Nz-2, Nz-2)
B = Array{ComplexF64, 2}(undef, Nz-2, Nz-2)

k = 2

T = zeros(Float64, (Nx, Nz))
T_out = zeros(Float64, (Nx, Nz))

T[:, k] .=  sin.(kx .* x)

# input = [fft(T[:, k], 1)[kx+1] for k in range(1, Nz)]
A[:, k-1] .= fft(T[:, 2:end-1], 1)[kx+1, :]


T_out .= implicit_solver(u, v, w, T)

# output = [fft(RHS[:, k], 1)[kx+1] for k in range(1, Nz)]
B[:, k-1] .= fft(T_out[:, 2:end-1], 1)[kx+1, :]
