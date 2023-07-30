function pressure_profile_for_params_unnormalised_x(xx::Float64, cfg::Config, params::Parameters)
    a1 = params.a1

    x_range = cfg.x_range

    β0 = maximum(
        [ 2 * a1 * magnetic_field_inner_iterate(xx, 0.0, cfg, params) for xx in x_range ]
    )

    pressure_profile(x) = β0 - 2 * a1 * magnetic_field_inner_iterate(x, 0.0, cfg, params)
                
    pressure_profile(xx)
end

function pressure_profile_for_params_unnormalised(cfg::Config, params::Parameters)
    x_range = cfg.x_range

    [ pressure_profile_for_params_unnormalised_x(xx, cfg, params) for xx in x_range ]
end

function pressure_profile_for_params_normalised_x(xx::Float64, cfg::Config, params::Parameters)
    p_vals = pressure_profile_for_params_unnormalised(cfg, params)
    M = maximum(abs.(p_vals))

    pressure_profile_for_params_unnormalised_x(xx, cfg, params) / M
end

function pressure_profile_for_params_normalised(cfg::Config, params::Parameters)

    p_vals = pressure_profile_for_params_unnormalised(cfg, params)

    M = maximum(abs.(p_vals))

    p_vals / M
end