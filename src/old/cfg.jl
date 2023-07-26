"""
Settings for the simulation
"""
mutable struct Configuration 
    α::Float64          # yep, see wang
    a1::Float64         # Another wang
    a2::Float64         # Another wang
    In1::Float64        # Current on the border
    In0::Tuple{Float64, Float64}        # Internal current (r0, In0)
    βv::Float64         # Pressure average
    x0::Float64         # Centre of Tokamak cross-section
    ρ::Float64          # Cross-section radius
    k::Float64          # Elongation of cross-section
    h::Int              # Simulation precision
    N::Int              # Eigenvalue precision
    x_range::StepRangeLen  # Values to compute for
    z_range::StepRangeLen  # Values to compute for
end

"""
Constructor for above. This should be used as opposed to
    manually creating a Configuration
"""
function Configuration(α::Float64, a1::Float64, a2::Float64, ρ::Float64, x0::Float64, h::Int)
    x_range = range(x0 - ρ, x0 + ρ, length=h)
    z_range = range(-ρ, ρ, length=h)

    Configuration(
        α,
        a1, 
        a2, 
        0.0, # void
        (0.0, 0.0), # void
        0.0, # void
        x0,
        ρ, 
        ρ, 
        h, 
        10, 
        x_range, 
        z_range
    )
end

function Configuration(In1::Float64, In0::Tuple{Float64, Float64}, βv::Float64,
    ρ::Float64, x0::Float64, h::Int, N::Int)
    x_range = range(x0 - ρ, x0 + ρ, length=h)
    z_range = range(-ρ, ρ, length=h)

    Configuration(
        0.0, # void
        0.0, # void
        0.0, # void
        In1,
        In0,
        βv,
        x0,
        ρ, 
        ρ, 
        h, 
        N, 
        x_range, 
        z_range
    )
end

"""
Calculated values for the given configuration.
"""
mutable struct Parameters
    aus::Vector{Float64}    
    ads::Vector{Float64}
    μns::Vector{Float64}    # 
    cns::Vector{Float64}    # 
    vls::Vector{Float64}    # 
    λnls::Matrix{Float64}   # 
end