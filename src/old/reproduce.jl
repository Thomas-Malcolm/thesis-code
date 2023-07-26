using Plots
using BenchmarkTools
using NLsolve

include("./computations.jl")
include("./current_profile.jl")
include("./cfg.jl")



"""
Calculate all values necessary for magnetic field
    and current calculations.
"""
function base_params(cfg::Configuration)::Parameters

    # Eigenvalues
    μns = generate_μn(
        cfg.N, cfg.x0
    )

    cns = generate_cn(
        cfg.x0, μns
    )

    vls = generate_vl(
        cfg.N, cfg.k
    )

    λnls = generate_eigenvalues(
        vls, μns
    )

    return Parameters(
        Vector([]), # void
        Vector([]), # void
        μns,
        cns,
        vls,
        λnls
    )
end

"""
Returns (a1, a2, α). Solves for them given current and pressure information.
"""
function solve_for_config(cfg, params)

    #=function f!(F, x)
        #=
        x[1] : a1
        x[2] : a2
        x[3] : α
        =#
        r0 = cfg.In0[1]
        In0 = cfg.In0[2]
        In1 = cfg.In1
        βv = cfg.βv
        F[1] = -In1 + current_for_params(x[1], x[2], x[3], 1, cfg, params)
        F[2] = -In0 + current_for_params(x[1], x[2], x[3], r0, cfg, params)
        F[3] = -βv + pressure_average(x[1], x[2], x[3], cfg, params)
    end=#

    function f_test(x)
        #=
        x[1] : a1
        x[2] : a2
        x[3] : α
        =#
        r0 = cfg.In0[1]
        In0 = cfg.In0[2]
        In1 = cfg.In1
        βv = cfg.βv
        f1 = -In1 + current_for_params(x[1], x[2], x[3], 1, cfg, params)
        f2 = -In0 + current_for_params(x[1], x[2], x[3], r0, cfg, params)
        f3 = -βv + pressure_average(x[1], x[2], x[3], cfg, params)
        return (f1, f2, f3)
    end

    f1, f2, f3 = f_test(([-0.04531; -1.0808; 2.1683]))
    println(f1)
    println(f2)
    println(f3)
    exit(1)

    # Don't have Jacobian (for now at least, is possible)
    r = nlsolve(f!, [-0.04531; -1.0808; 2.1683])

    println(r)

    if !converged(r)
        println("Failed to converge for a1, a2, α")
        exit(1)
    end

    return r.zero[1], r.zero[2], r.zero[3]
end

function run_simulation(In1::Float64, In0::Tuple{Float64, Float64}, βv::Float64, ρ::Float64, x0::Float64, h::Int, N::Int)

    cfg = Configuration(
        In1,
        In0,
        βv,
        ρ,
        x0,
        h,
        N
    )

    params = base_params(cfg)

    #println("Solving for a1, a2 and α")
    #a1, a2, α = solve_for_config(cfg, params)
    
    #cfg.a1 = a1
    #cfg.a2 = a2
    #cfg.α = α

    cfg.a1 = -0.04531
    cfg.a2 = -1.0808
    cfg.α = 2.1683

    params.aus = aun.(1:N, Ref(cfg), Ref(params.μns), Ref(params.cns))
    params.ads = adn.(1:N, Ref(cfg), Ref(params.μns), Ref(params.cns))

    println("Generating magnetic field...")
    mag_field = @time magnetic_field(cfg, params)

    println("Calculating current profile...")
    ρ_vals = range(-cfg.ρ, cfg.ρ, length=cfg.h)
    current_values = @time current_profile.(ρ_vals, Ref(cfg), Ref(params))

    cont = contour(cfg.x_range, cfg.z_range, mag_field)
    plot!(cont, cfg.x_range, current_values)
    savefig(cont, "banana_boat.pdf")

end

function plot_current(
    a1::Float64,
    a2::Float64,
    α::Float64,
    ρ::Float64,
    x0::Float64,
    h::Int,
    N::Int
)

    cfg = Configuration(
        α,
        a1,
        a2,
        ρ,
        x0,
        h
    )

    params = base_params(cfg)
    params.aus = aun.(1:N, Ref(cfg), Ref(params.μns), Ref(params.cns))
    params.ads = adn.(1:N, Ref(cfg), Ref(params.μns), Ref(params.cns))

    println("Calculating current profile...")
    ρ_vals = range(-cfg.ρ, cfg.ρ, length=cfg.h)
    current_values = @time current_profile.(ρ_vals, Ref(cfg), Ref(params))

    p = plot(ρ_vals, current_values)
    savefig(p, "banana_boat.pdf")
end 

#=
plot_current(
    2.1683,
    -0.04531,
    -1.0808,
    1.0,
    5.0,
    200,
    10
)
=#


run_simulation(
    0.33,
    (0.7, 0.5),
    0.01,
    1.0,
    5.0,
    400,
    50
)

#=

function run_simulation()
    cfg = Configuration(
        2.1683,
        -0.04531,
        -1.0808,
        1.0,
        5.0,
        500
    )
    #=cfg = Configuration(
        5.566,
        0.01,
        3.1,
        1.0,
        5.0,
        2000
    )=#
    println(cfg)

    println("Producing parameters...")
    params = produce_parameters(cfg)

    println("Generating magnetic field...")
    mag_field = @time magnetic_field(cfg, params)

    println("Calculating current profile...")
    ρ_vals = range(-cfg.ρ, cfg.ρ, length=cfg.h)
    current_values = @time current_profile.(ρ_vals, Ref(cfg), Ref(params))

    cont = contour(cfg.x_range, cfg.z_range, mag_field)
    plot!(cont, cfg.x_range, current_values)
    savefig(cont, "figure2-c.pdf")
end
=#

