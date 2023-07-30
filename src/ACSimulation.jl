module ACSimulation

using Bessels
using Roots

export Config
export Parameters
export Datum
export Data

export initial_calculations
export auns_for_params

export magnetic_field_for_params
export magnetic_field_inner_iterate
export ψ

export current_profile_for_params_unnormalised_x
export current_profile_for_params_unnormalised
export current_profile_for_params_normalised_x
export current_profile_for_params_normalised

export pressure_profile_for_params_unnormalised_x
export pressure_profile_for_params_unnormalised
export pressure_profile_for_params_normalised_x
export pressure_profile_for_params_normalised

export CurrentDa1
export CurrentDa2
export CurrentDalpha

export residual
export residuals

export residualDa1
export residualDa2
export residualDalpha


mutable struct Config
    x0::Float64
    ρ::Float64
    k::Float64
    h::Int
    N::Int

    x_range
    z_range

    vls::Vector{Float64}
    μns::Vector{Float64}
    cns::Vector{Float64}
    λnls::Matrix{Float64}
    auns::Vector{Float64}
    adns::Vector{Float64}

    function Config(
        ;
        x0::Float64 = 5.0,
        ρ::Float64 = 1.0,
        k::Float64 = ρ,
        h::Int = 200,
        N::Int = 8,
    )

    new_cfg = new(
        x0,
        ρ,
        k,
        h,
        N,
        range(x0 - ρ, x0 + ρ, length = h),
        range(-ρ, ρ, length = h),

        Vector{Float64}(),
        Vector{Float64}(),
        Vector{Float64}(),
        Matrix{Float64}(undef, N, N),
        Vector{Float64}(),
        Vector{Float64}()
    )

    initial_calculations(new_cfg)
    new_cfg
    end
end

# Container for a1, a2 and α values
struct Parameters
    a1::Float64
    a2::Float64
    α::Float64
end

# Pairing of x val and its expected value
struct Datum
    x::Float64
    val::Float64
end

struct Data
    d::Vector{Datum}

    function Data(
        arr::Vector{Tuple{Float64, Float64}}
    )

        d_vec = [
            Datum(x, val) for (x, val) in arr
        ]

        new(
            d_vec
        )
    end
end

function initial_calculations(cfg::Config)

    h = cfg.h
    N = cfg.N
    k = cfg.k

    x0 = cfg.x0
    ρ = cfg.ρ

    x_range = cfg.x_range
    z_range = cfg.z_range

    # Vls
    cfg.vls = [
        (l + 0.5) * π / k for l in 0:N-1
    ]
    vls = cfg.vls

    # μns
    zero_function(μn) = -besselj1(μn * (x0 + 1)) * ( bessely1(μn * (x0 - 1)) / besselj1(μn * (x0 - 1)) ) + bessely1(μn * (x0 + 1))

    zeros = find_zeros(
        zero_function, (0, 2π*N)
    )

    true_zeros = findall(
        μn -> abs(zero_function(μn)) < 10e-10, zeros
    )

    cfg.μns = zeros[splice!(true_zeros, 1:N)]
    μns = cfg.μns

    # cns
    cn(μn) = -bessely1(μn * (x0 - 1)) / besselj1(μn * (x0 - 1))

    cfg.cns = cn.(μns)
    cns = cfg.cns

    # λnls
    cfg.λnls = [ vls[l]^2 + μns[n]^2 for n in 1:N, l in 1:N ]
    λnls = cfg.λnls

    # adns
    edn(x, n) = -0.5 * x^2 * bessely0(μns[n] * x) * bessely(2, μns[n] * x) - 0.5 * cns[n]^2 * x^2 * besselj0(μns[n] * x) * besselj(2, μns[n] * x) + cns[n] * (
        0.5 * x^2 * (
            besselj0(μns[n] * x) - besselj(2, μns[n] * x)
        )
        * bessely0(μns[n] * x)
        -
        (1/μns[n]) * x * besselj0(μns[n] * x) * bessely1(μns[n] * x)
    )

    cfg.adns = [ edn(x0 + 1, n) - edn(x0 - 1, n) for n in 1:N ]
end

"""
Given a specific a1, a2, α value, generate the au(n) values.
"""
function auns_for_params(cfg::Config, params::Parameters)
    
    a1 = params.a1
    a2 = params.a2
    α = params.α
    
    μns = cfg.μns
    cns = cfg.cns
    x0 = cfg.x0
    N = cfg.N
    
    eun(x, n) = (1 / μns[n]) * ( 
        a1 * x^2 * (
            cns[n] * besselj(2, μns[n] * x) + bessely(2, μns[n] * x)
        )
        +
        a2 * (
            cns[n] * besselj0(μns[n] * x) + bessely0(μns[n] * x)
        )
    )

    auns = [ (eun(x0 + 1, n) - eun(x0 - 1, n)) for n in 1:N ]
    auns
end

function magnetic_field_inner_iterate(xx::Float64, zz::Float64, cfg::Config, params::Parameters)
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

    inner_iterate(xx, zz)
end

ψ(xx::Float64, zz::Float64, cfg::Config, params::Parameters) = magnetic_field_inner_iterate(xx, zz, cfg, params)

function magnetic_field_for_params(cfg::Config, params::Parameters)
    x_range = cfg.x_range
    z_range = cfg.z_range

    mag_field = [
        magnetic_field_inner_iterate(x,z, cfg, params) for z in z_range, x in x_range
    ]

    mag_field
end

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

## Derivatives

### Current Derivatives

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
    k = cfg.k
    N = cfg.N

    tmpVal(n) = (1/μns[n]) * (
        (x0+1)^2 * (
            cns[n] * besselj(2, μns[n] * (x0 + 1)) + bessely(2, μns[n] * (x0 + 1))
        ) 
        -
        (x0-1)^2 * (
            cns[n] * besselj(2, μns[n] * (x0 - 1)) + bessely(2, μns[n] * (x0 - 1))
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
    k = cfg.k
    N = cfg.N

    tmpVal(n) = (1/μns[n]) * (
        (
            cns[n] * besselj0(μns[n] * (x0 + 1)) + bessely0(μns[n] * (x0 + 1))
        ) 
        -
        (
            cns[n] * besselj0(μns[n] * (x0 - 1)) + bessely0(μns[n] * (x0 - 1))
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

# Data Based Functions

"""
Calculates the residual value from the expected (`d`), and our calculated
    based on the Config and provided `a1, a2, α` parameters
"""
function residual(cfg::Config, params::Parameters, d::Datum)
    M = maximum_current_for_params(cfg, params)

    d.val - (1/M) * current_profile_for_params_unnormalised_x(d.x, cfg, params)
end

"""
Calculates the total residual squared value from a set of data points. 
"""
function residuals(cfg::Config, params::Parameters, ds::Data)

    sum(residual(cfg, params, d) for d in ds.d)
end

## Residual Derivatives

"""
Partial derivative of the `residuals` function with respect to a1, for a given 
    set of data
"""
function residualDa1(cfg::Config, params::Parameters, ds::Data)
    
    M = maximum_current_for_params(cfg, params)

    sum(
        (-2 / M) * (
            d.val - (1/M) * current_profile_for_params_unnormalised_x(d.x, cfg, params)
        ) * CurrentDa1(d.x, cfg, params)
        for d in ds.d
    )
end

"""
Partial derivative of the `residuals` function with respect to a2, for a given 
    set of data
"""
function residualDa2(cfg::Config, params::Parameters, ds::Data)
    M = maximum_current_for_params(cfg, params)

    sum(
        (-2 / M) * (
            d.val - (1/M) * current_profile_for_params_unnormalised_x(d.x, cfg, params)
        ) * CurrentDa2(d.x, cfg, params)
        for d in ds.d
    )
end

"""
Partial derivative of the `residuals` function with respect to α, for a given 
    set of data
"""
function residualDalpha(cfg::Config, params::Parameters, ds::Data)
    M = maximum_current_for_params(cfg, params)

    sum(
        (-2 / M) * (
            d.val - (1/M) * current_profile_for_params_unnormalised_x(d.x, cfg, params)
        ) * CurrentDalpha(d.x, cfg, params)
        for d in ds.d
    )
end



end # Module end