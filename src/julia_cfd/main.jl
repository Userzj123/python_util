module adjoint_lesgo

include("meshgrid/mesh.jl")
include("params.jl")
include("scalar/scalar.jl")
include("momentum/velocity.jl")
using .mesh
using .params
using .scalar_solver
using .velocity_solver


struct vars
    velocity :: velocity_vars
    scalar :: scalar_vars
    conf :: configs 
    grid :: meshgrid 
end



end 