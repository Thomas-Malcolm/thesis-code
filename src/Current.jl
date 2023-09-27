

function current_profile_for_params_unnormalised_x(xx::Float64, cfg::Config, params::Parameters)
    a1 = params.a1
    a2 = params.a2
    α = params.α
    
    
    auns = auns_for_params(cfg, params)
    adns = cfg.adns
    μns = cfg.μns
    vls = cfg.vls
    λnls = cfg.λnls
    cns = cfg.cns
    k = cfg.k
    N = cfg.N
    
    inner_iterate(x, z) = x * sum(
        sum(
            (((-1)^(l-1) * 2 * auns[n]) / (k * vls[l] * adns[n] * (α^2 - λnls[n,l])))
            *
            (cns[n] * besselj1(μns[n] * x) + bessely1(μns[n] * x))
            *
            cos(vls[l] * z)
            for l in 1:N
        ) for n in 1:N
    )
    
    current_profile(x) = -a1 * x + a2/x + (1/x)*(α^2) * inner_iterate(x,0)
    current_profile(xx)
end

function current_profile_for_params_unnormalised(cfg::Config, params::Parameters)
    x_range = cfg.x_range

    [ current_profile_for_params_unnormalised_x(xx, cfg, params) for xx in x_range ]
end

function current_profile_for_params_normalised_x(xx::Float64, cfg::Config, params::Parameters)
    M = maximum_current_for_params(cfg, params)

    current_profile_for_params_unnormalised_x(xx, cfg, params) / M
end

function current_profile_for_params_normalised(cfg::Config, params::Parameters)
    curr_vals = current_profile_for_params_unnormalised(cfg, params)

    M = maximum_current_for_params(cfg, params)

    curr_vals / M
end

function maximum_current_for_params(cfg::Config, params::Parameters)
    curr_vals = current_profile_for_params_unnormalised(cfg, params)

    maximum(abs.(curr_vals))
end

function toroidal_current_density_profile_x(xx::Float64, cfg::Config, params::Parameters)
    jϕ = current_profile_for_params_unnormalised_x(xx, cfg, params)

    μ0 = cfg.μ0
    B0 = cfg.B0
    a = cfg.ρ

    (jϕ * B0) / (μ0 * a)
end

function toroidal_current_density_profile(cfg::Config, params::Parameters)
    x_range = cfg.x_range

    [ toroidal_current_density_profile_x(xx, cfg, params) for xx in x_range ]
end


# Current Derivatives

"""
Partial derivative of unnormalised current with respect to a1
"""
function CurrentDa1(xx::Float64, cfg::Config, params::Parameters)

    α = params.α

    μns = cfg.μns
    cns = cfg.cns
    vls = cfg.vls
    λnls = cfg.λnls
    adns = cfg.adns

    x0 = cfg.x0
    ρ = cfg.ρ
    k = cfg.k
    N = cfg.N

    tmpVal(n) = (1/μns[n]) * (
        (x0+1)^2 * (
            cns[n] * besselj(2, μns[n] * (x0 + ρ)) + bessely(2, μns[n] * (x0 + ρ))
        ) 
        -
        (x0-1)^2 * (
            cns[n] * besselj(2, μns[n] * (x0 - ρ)) + bessely(2, μns[n] * (x0 - ρ))
        )
    )

    -xx + α^2 * sum(
        sum(
            (
                ( (-1)^(l-1) * 2 * tmpVal(n)) 
                / 
                ( k * vls[l] * adns[n] * (α^2 - λnls[n,l]) )
            )
            *
            (cns[n] * besselj1(μns[n] * xx) + bessely1(μns[n] * xx) )
            for l in 1:N
        ) for n in 1:N
    )
end

"""
Partial derivative of unnormalised current with respect to a2
"""
function CurrentDa2(xx::Float64, cfg::Config, params::Parameters)

    α = params.α

    μns = cfg.μns
    cns = cfg.cns
    vls = cfg.vls
    λnls = cfg.λnls
    adns = cfg.adns

    x0 = cfg.x0
    ρ = cfg.ρ
    k = cfg.k
    N = cfg.N

    tmpVal(n) = (1/μns[n]) * (
        (
            cns[n] * besselj0(μns[n] * (x0 + ρ)) + bessely0(μns[n] * (x0 + ρ))
        ) 
        -
        (
            cns[n] * besselj0(μns[n] * (x0 - ρ)) + bessely0(μns[n] * (x0 - ρ))
        )
    )

    (1/xx) + α^2 * sum(
        sum(
            (
                ( (-1)^(l-1) * 2 * tmpVal(n)) 
                / 
                ( k * vls[l] * adns[n] * (α^2 - λnls[n,l]) )
            )
            *
            (cns[n] * besselj1(μns[n] * xx) + bessely1(μns[n] * xx) )
            for l in 1:N
        ) for n in 1:N
    )
end

"""
Partial derivative of unnormalised current with respect to α
"""
function CurrentDalpha(xx::Float64, cfg::Config, params::Parameters)

    a1 = params.a1
    a2 = params.a2
    α = params.α

    μns = cfg.μns
    cns = cfg.cns
    vls = cfg.vls
    λnls = cfg.λnls
    adns = cfg.adns
    auns = auns_for_params(cfg, params)

    N = cfg.N
    k = cfg.k

    part1 = 2 * α * (1/xx) * ψ(xx, 0.0, cfg, params)

    part2 = α^2 * sum(
        sum(
            (
                ((-1)^(l) * 4 * auns[n] * α)
                / 
                (k * vls[l] * adns[n] * (α^2 - λnls[n,l])^2)
            )
            *
            (cns[n] * besselj1(μns[n] * xx) + bessely1(μns[n] * xx))
            for l in 1:N
        ) for n in 1:N
    )

    part1 + part2    
end