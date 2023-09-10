module interp

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


end