
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

function magnetic_field_for_params_and_range(cfg::Config, params::Parameters, x_range, z_range)
    mag_field = [
        magnetic_field_inner_iterate(x, z, cfg, params) for z in z_range, x in x_range
    ]

    mag_field
end
