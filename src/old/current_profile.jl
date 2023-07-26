include("cfg.jl")
include("computations.jl")

"""
Equation 15c from Wang
"""
function current_profile(ρ:: Float64, cfg::Configuration, params::Parameters)

    x0 = cfg.x0
    k = cfg.k
    a1 = cfg.a1
    a2 = cfg.a2
    α = cfg.α
    N = cfg.N
    aus = params.aus
    ads = params.ads
    vls = params.vls
    cns = params.cns
    λnls = params.λnls
    μns = params.μns
    
    x0p = x0 + ρ
    x0n = x0 - ρ

    sumOne = (
        -(π * x0) / (8 * ρ^2 * (1 + k)^2) * k * ρ 
        * ( 
            (a1 * x0p^2 - 2 * a2 * log(x0p)) 
            - 
            (a1 * x0n^2 - 2 * a2 * log(x0n)) 
        )
    )

    sumTwo = (
        -(π * x0) / (8 * ρ^2 * (1 + k)^2) * α^2
        * sum(
            sum(
                (
                    (-1)^(l-1) * 4 * aus[n] * sin(vls[l] * k * ρ)
                )
                /
                (
                    k * μns[n] * vls[l]^2 + ads[n] * (α^2 - λnls[n,l])
                )
                *
                (
                    (
                        cns[n] * besselj0(μns[n] * x0p) + bessely0(μns[n] * x0p)
                    )
                    -
                    (
                        cns[n] * besselj0(μns[n] * x0n) + bessely0(μns[n] * x0n)
                    )
                ) for l in 1:N
            ) for n in 1:N
        ) 
    )

    return sumOne + sumTwo
end

function current_for_params(a1, a2, α, ρ, cfg, params)
    x0 = cfg.x0
    k = cfg.k
    N = cfg.N
    vls = params.vls
    cns = params.cns
    λnls = params.λnls
    μns = params.μns
    ads = adn.(1:N, Ref(cfg), Ref(μns), Ref(cns))
    
    x0p = x0 + ρ
    x0n = x0 - ρ

    sumOne = (
        -(π * x0 * k * ρ ) / (8 * ρ^2 * (1 + k)^2) 
        * ( 
            (a1 * x0p^2 - 2 * a2 * log(x0p))
            - 
            (a1 * x0n^2 - 2 * a2 * log(x0n)) 
        )
    )

    sumTwo = (
        -(π * x0 * α^2) / (8 * ρ^2 * (1 + k)^2) 
        * sum(
            sum(
                (
                    (-1)^(l-1) * 4 * aun(n, a1, a2, μns, cns, cfg) * sin(vls[l] * k * ρ)
                )
                /
                (
                    k * μns[n] * vls[l]^2 + ads[n] * (α^2 - λnls[n,l])
                )
                *
                (
                    (
                        cns[n] * besselj0(μns[n] * x0p) + bessely0(μns[n] * x0p)
                    )
                    -
                    (
                        cns[n] * besselj0(μns[n] * x0n) + bessely0(μns[n] * x0n)
                    )
                ) for l in 1:N
            ) for n in 1:N
        ) 
    )

    return sumOne + sumTwo
end

function pressure_average(a1, a2, α, cfg, params)
    x0 = cfg.x0
    k = cfg.k
    N = cfg.N
    h = cfg.h
    ρ = cfg.ρ
    vls = params.vls
    cns = params.cns
    λnls = params.λnls
    μns = params.μns
    ads = adn.(1:N, Ref(cfg), Ref(μns), Ref(cns))
    
    x0p = x0 + ρ
    x0n = x0 - ρ

    magnetic_lines = magnetic_field(
        α,
        cfg,
        params
    )

    β0 = 2 * a1 * minimum(
        minimum(
            magnetic_lines[i,:] for i in 1:h
        )
    )

    pressure_val = (
        (a1 / (2 * k * x0) )
        *
        sum(
            sum(
                (
                    4 * aun(n, a1, a2, μns, cns, cfg)
                )
                /
                (
                    k * μns[n] * vls[l]^2 + ads[n] * (α^2 - λnls[n,l])
                )
                *
                (
                    (
                        cns[n] * x0p^2 * besselj(2, μns[n] * x0p) + bessely(2, μns[n] * x0p)
                    )
                    -
                    (
                        cns[n] * x0n^2 * besselj(2, μns[n] * x0n) + bessely(2, μns[n] * x0n)
                    )
                ) for l in 1:N
            ) for n in 1:N
        ) 
    )

    return β0 - pressure_val
end