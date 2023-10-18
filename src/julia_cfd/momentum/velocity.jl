module velocity_solver

export velocity_vars, interp_veloctiy_vars


struct velocity_vars
    u :: vartype
    v :: vartype
    w :: vartype
end


struct interp_veloctiy_vars
    u_iijk :: vartype
    v_ijjk :: vartype
    w_ijkk :: vartype
end




end