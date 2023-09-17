module main
using LinearAlgebra
using FFTW

export explicit_solver





    function explicit_solver(u, v, w, T)
        (Nx, Nz) = size(u)
        adv = zeros(Float64, (Nx, Nz))
        adv[:, 2:end-1] .= u[:, 2:end-1] .* (circshift(T[:, 2:end-1],[-1, 0]) .- circshift(T[:, 2:end-1],[1, 0]))./dx./2

        diff = zeros(Float64, (Nx, Nz))
        diff[:, 2:end-1] .= nu .*( circshift(T[:, 2:end-1],[-1,0]) .- 2 .*T[:, 2:end-1] .+ circshift(T[:, 2:end-1],[1,0]))./dx./dx .+ 
        nu .*(T[:, 3:end] .- 2 .*T[:, 2:end-1] .+ T[:, 1:end-2] ) ./dz./dz

        T_out = zeros(Float64, (Nx, Nz))
        T_out .= T .- dt .* adv .+ dt .* diff

        return T_out
    end


end