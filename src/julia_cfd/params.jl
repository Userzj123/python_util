module params

using TOML

export varprec, vartype, boundaries, configs, read_conf

const varprec = Float64
const vartype = Union{Array{varprec, 3}, Array{varprec, 4}}


struct boundary
    type :: Int64       # 1 - Neumann  0 - Dirichlet
    value :: Float64    # if type = 1 -- scalar flux, type = 0 -- prescribed scalar
end

struct boundaries
    lbc :: boundary
    ubc :: boundary
end


struct configs
    domain :: Dict{String, Any}
    scalar :: Dict{String, Any}
    scalar_bcs :: boundaries
end

# struct configs
#     Nx :: Int64
#     Ny :: Int64
#     Nz :: Int64

#     Pr :: Float64
#     dt :: Float64
#     scalar_bc :: boundaries
# end


    function read_conf(fname::String)
        config::Dict = TOML.parsefile(fname)
        lbc = boundary(config["SCALAR"]["lbc_scalar"], config["SCALAR"]["lbc_value"])
        ubc = boundary(config["SCALAR"]["ubc_scalar"], config["SCALAR"]["ubc_value"])
        bcs = boundaries(lbc, ubc)
        return configs(config["DOMAIN"], config["SCALAR"], bcs)
    end


end