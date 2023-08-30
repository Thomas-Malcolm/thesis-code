


mutable struct Config
    x0::Float64                 # Centre of cross-section (m)
    ρ::Float64                  # Radius of cross-section (m)
    k::Float64                  # Elongation of cross-section (vertical radius) (m)
    h::Int                      # Mesh precision
    N::Int                      # Eigenvalue precision
    B0::Float64                 # On-axis magnetic field strength (T)   
    μ0::Float64                 # free-space permeability (N * A^(-2))

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
        B0::Float64 = 1.0,
        μ0::Float64 = 1.25663706212e-6
    )

    new_cfg = new(
        x0,
        ρ,
        k,
        h,
        N,
        B0,
        μ0,
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


function initial_calculations(cfg::Config)

    h = cfg.h
    N = cfg.N
    k = cfg.k

    # Radius information
    x0 = cfg.x0     # Major radius
    ρ = cfg.ρ       # Minor radius

    x_range = cfg.x_range
    z_range = cfg.z_range

    # Vls
    cfg.vls = [
        (l + 0.5) * π / k for l in 0:N-1
    ]
    vls = cfg.vls

    # μns
    zero_function(μn) = -besselj1(μn * (x0 + ρ)) * ( bessely1(μn * (x0 - ρ)) / besselj1(μn * (x0 - ρ)) ) + bessely1(μn * (x0 + ρ))

    zeros = find_zeros(
        zero_function, (0, 2π*N*10) # hope this is enough
    )

    true_zeros = findall(
        μn -> abs(zero_function(μn)) < 10e-8, zeros
    )

    cfg.μns = zeros[splice!(true_zeros, 1:N)]
    μns = cfg.μns

    # cns
    cn(μn) = -bessely1(μn * (x0 - ρ)) / besselj1(μn * (x0 - ρ))

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

    cfg.adns = [ edn(x0 + ρ, n) - edn(x0 - ρ, n) for n in 1:N ]
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
    ρ = cfg.ρ
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

    auns = [ (eun(x0 + ρ, n) - eun(x0 - ρ, n)) for n in 1:N ]
    auns
end
