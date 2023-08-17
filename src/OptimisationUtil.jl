
"""
Calculates the residual value from the expected (`d`), and our calculated
    based on the Config and provided `a1, a2, α` parameters
"""
function residual(cfg::Config, params::Parameters, d::Datum)
    M = maximum_current_for_params(cfg, params)

    d.val - (1/M) * current_profile_for_params_unnormalised_x(d.x, cfg, params)
end

function un_residual(cfg::Config, params::Parameters, d::Datum)
    M = maximum_current_for_params(cfg, params)

    d.val - current_profile_for_params_unnormalised_x(d.x, cfg, params)
end

"""
Calculates the total residual squared value from a set of data points. 
"""
function residuals(cfg::Config, params::Parameters, ds::Data)

    sum(residual(cfg, params, d)^2 for d in ds.d)
end

function un_residuals(cfg::Config, params::Parameters, ds::Data)

    sum(un_residual(cfg, params, d)^2 for d in ds.d)
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

function residualUnDa1(cfg::Config, params::Parameters, ds::Data)

    sum(
        -2 * (
            d.val - current_profile_for_params_unnormalised_x(d.x, cfg, params)
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

function residualUnDa2(cfg::Config, params::Parameters, ds::Data)

    sum(
        -2 * (
            d.val - current_profile_for_params_unnormalised_x(d.x, cfg, params)
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

function residualUnDalpha(cfg::Config, params::Parameters, ds::Data)

    sum(
        -2 * (
            d.val - current_profile_for_params_unnormalised_x(d.x, cfg, params)
        ) * CurrentDalpha(d.x, cfg, params)
        for d in ds.d
    )
end
