module scalar

include("./mesh.jl")
include("./params.jl")
using .mesh
using .params

export scalar_forward!


struct scalar_vars
    theta :: vartype
    source :: vartype
    RHS :: vartype
    advection :: vartype
    diffusion :: vartype
end


struct velocity_vars
    u :: vartype
    v :: vartype
    w :: vartype
end

struct vars
    velocity :: velocity_vars
    scalar :: scalar_vars
    conf :: configs 
    grid :: meshgrid 
end

struct interp_veloctiy
    u_iijk :: vartype
    v_ijjk :: vartype
    w_ijkk :: vartype
end

struct interp_scalar
    theta_iijk :: vartype
    theta_ijjk :: vartype
    theta_ijkk :: vartype
end

struct interp_vars
    velocity :: interp_veloctiy
    scalar :: interp_scalar
end

    function scalar_forward!(input::vars, )

        interp_result::interp_vars = interpolation(input)

        advection::Array = get_adv!(input, interp_result)
        # diffusion::Array = get_diff!(input)

        # RHS:: Array = - advection + diffusion + source

        # # Time scalar_forward
        # c = c +  dt * RHS
    end


    function get_adv(input::vars, interp_result::interp_vars)
        # Working On Sep 9

        advection = Array{Float64, 3}(undef, input.conf.domain["nx"], input.conf.domain["ny"], input.conf.domain["nz"])
        for i in Range(Nx)
            for j in Range(Ny)
                for k in Range(Nz)
                    advection[i, j, k] = (
                        interp_result.velocity.u_iijk[i+1, j, k]*input.grid.S.Sx[i, j, k]*interp_result.scalar.theta_iijk[i+1, j, k]    - 
                        interp_result.velocity.u_iijk(ji, jj, jk)*input.grid.S.Sx(jk)*interp_result.scalar.theta_iijk(i, j, k)        + 
                        interp_result.velocity.v_ijjk(ji, jj+1, jk)*input.grid.S.Sy(jk)*interp_result.scalar.theta_ijjk(ji, jj+1, jk)    - 
                        interp_result.velocity.v_ijjk(ji, jj, jk)*input.grid.S.Sy(jk)*interp_result.scalar.theta_ijjk(ji, jj, jk)        + 
                        interp_result.velocity.w_ijkk(ji, jj, jk+1)*input.grid.S.Sz(jk)*interp_result.scalar.theta_ijkk(ji, jj, jk+1)    - 
                        interp_result.velocity.w_ijkk(ji, jj, jk)*input.grid.S.Sz(jk)*interp_result.scalar.theta_ijkk(ji, jj, jk)
                    )/input.grid.V(jk)
                end 
            end
        end
        return advection
    end

    function get_diff(input::vars)
    end


    function interpolation(input::vars)::interp_vars
        u_iijk = interp_ijk2iijk(input.velocity.u)
        v_ijjk = interp_ijk2ijjk(input.velocity.v)
        w_ijkk = input.velocity.w

        vel = interp_veloctiy(u_iijk, v_ijjk, w_ijkk)

        theta_iijk = interp_ijk2iijk(input.scalar.theta)
        theta_ijjk = interp_ijk2ijjk(input.scalar.theta)
        theta_ijkk = interp_ijk2ijkk(input.scalar.theta, input.conf.scalar_bc, input.grid.ijk, input.grid.ijkk)
        c = interp_scalar(theta_iijk, theta_ijjk, theta_ijkk)
        return interp_vars(vel, c)
    end


    function interp_ijk2iijk(f::vartype)::vartype
        Nx, Ny, Nz = size(f)
        f_ii = Array{varprec, 3}(undef, Nx+1, Ny, Nz)
        for jx in 1:Nx+1
            nhalf = mod(jx - 1, Nx)
            phalf = mod(jx, Nx)
            if (nhalf==0) nhalf = Nx end # Periodic in x direction
            if (phalf==0) phalf = Nx end

            f_ii[jx, :, :] .= (f[phalf, :, :] .+ f[nhalf, :, :])./2
        end
        return f_ii
    end


    function interp_ijk2ijjk(f::vartype)::vartype
        Nx, Ny, Nz = size(f)
        f_jj = Array{varprec, 3}(undef, Nx, Ny+1, Nz)
        for jy in 1:Ny+1 
            nhalf = mod(jy - 1, Ny)
            phalf = mod(jy, Ny)
            if (nhalf==0) nhalf = Ny end # Periodic in y direction
            if (phalf==0) phalf = Ny end

            f_jj[:, jy, :] .= (f[:, phalf, :] .+ f[:, nhalf, :])./2
        end
        return f_jj
    end

    function interp_ijk2ijkk(f::vartype, bcs::boundaries, ijk::coordinates, ijkk::coordinates)::vartype
        
        zc = ijk.z
        z = ijkk.z
        Nx, Ny, Nz = size(f)
        f_kk = Array{varprec, 3}(undef, Nx, Ny, Nz+1)
        for jz in 2:Nz 
            f_jj[:, :, jz] .= (f[:, :, jz] .- f[:, :, jz-1]) ./ 
                (zc[jz] - zc[jz-1]) .* (z[jz] - zc[jz-1]) .+ f[:, :, jz-1]
        end

        if bcs.lbc.type == 1
            f_jj[:, :, 1] .= f[:, :, 1] .- 0.5 * bcs.lbc.value*zc[1]*2
        elseif bcs.lbc.type == 0
            f_jj[:, :, 1] .= bcs.lbc.value
        end 

        if bcs.ubc.type == 1
            f_jj[:, :, Nz+1] .= f[:, :, Nz] .+ 0.5 * bcs.ubc.value * (z[Nz+1] - zc[Nz]) * 2
        elseif bcs.ubc.type == 0
            f_jj[:, :, Nz+1] .= bcs.ubc.value
        end

        return f_jj
    end



end