module scalar





struct coordinates(x, y, z)
    x :: Vector
    y :: Vector
    z :: Vector
    Nx :: Int = length(x)
    Ny :: Int = length(y)
    Nz :: Int = length(z)
end

struct mesh
    xx :: Array
    yy :: Array
    zz :: Array
end

global u::Array, v::Array, w::Array

global c_RHS::Array, c_advection::Array, c_diffusion::Array

global c::Array


function init(domain, dims)

end


function udpate_adv()

    for i in Range(Nx)
        for j in Range(Ny)
            for k in Range(Nz)
        c_advections = (u_ihalf(ji+1, jj, jk)*Sx(jk)*theta_ihalf(ji+1, jj, jk)  &
        - u_ihalf(ji, jj, jk)*Sx(jk)*theta_ihalf(ji, jj, jk)                    &
        + v_jhalf(ji, jj+1, jk)*Sy(jk)*theta_jhalf(ji, jj+1, jk)                &
        - v_jhalf(ji, jj, jk)*Sy(jk)*theta_jhalf(ji, jj, jk)                    &
        + w(ji, jj, jk+1)*Sz(jk)*theta_kw(ji, jj, jk+1)                           &
        - w(ji, jj, jk)*Sz(jk)*theta_kw(ji, jj, jk))/Volume(jk)
            end 
        end
    end

    # for i, x in enumerate(ldata.coords[0]):
    #     for j, y in enumerate(ldata.coords[1]):
    #         for k, z in enumerate(ldata.kw_coords[2][:-1]):
    #             tmp[i, j, k] = (
    #                 (ldata.data['u_ihalf'][i+1, j, k] * ldata.data['theta_ihalf'][nk, i+1, j, k]
    #                     -
    #                 ldata.data['u_ihalf'][i, j, k] * ldata.data['theta_ihalf'][nk, i, j, k]
    #                  )/(ldata.ihalf_coords[0][i+1] - ldata.ihalf_coords[0][i])
    #                 +
    #                 (ldata.data['v_jhalf'][i, j+1, k] * ldata.data['theta_jhalf'][nk, i, j+1, k]
    #                     -
    #                 ldata.data['v_jhalf'][i, j, k] * ldata.data['theta_jhalf'][nk, i, j, k]
    #                  )/(ldata.jhalf_coords[1][j+1] - ldata.jhalf_coords[1][j])
    #                 +
    #                 (ldata.data['w_kw'][i, j, k+1] * ldata.data['theta_kw'][nk, i, j, k+1]
    #                     -
    #                 ldata.data['w_kw'][i, j, k] * ldata.data['theta_kw'][nk, i, j, k]
    #                  )/(dz * gradient_str_func(Lz = ldata.domain[2], z = ldata.coords[2], str_factor=1.5)[0])
    #             ) 
end

end