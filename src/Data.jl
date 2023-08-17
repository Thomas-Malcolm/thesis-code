

# Pairing of x val and its expected value
struct Datum
    x::Float64
    val::Float64
end

function -(d1::Datum, d2::Datum)
    @assert(d1.x == d2.x)

    Datum(d1.x, d1.val - d2.val)
end

function -(d::Datum, val::Float64)
    Datum(d.x, d.val - val)
end

function +(d1::Datum, d2::Datum)
    @assert(d1.x == d2.x)

    Datum(d1.x, d1.val + d2.val)
end

function +(d::Datum, val::Float64)
    Datum(d.x, d.val + val)
end

function /(d::Datum, denom::Float64)
    @assert(denom != 0.0)

    Datum(d.x, d.val / denom)
end

function /(d::Datum, denom::Int64)
    d / float(denom)
end

function *(d::Datum, denom::Float64)
    @assert(denom != 0.0)

    Datum(d.x, d.val * denom)
end

function *(d::Datum, denom::Int64)
    d * float(denom)
end

function -(d1v::Vector{Datum}, d2v::Vector{Datum})

    Vector{Datum}(
        [
            d1 - d2v[j]
            for (j, d1) in pairs(d1v)
        ]
    )
end

function -(d1v::Vector{Datum}, val::Float64)

    Vector{Datum}(
        [
            d1 - val
            for d1 in d1v
        ]
    )
end

function +(d1v::Vector{Datum}, d2v::Vector{Datum})

    Vector{Datum}(
        [
            d1 + d2v[j]
            for (j, d1) in pairs(d1v)
        ]
    )
end

function +(d1v::Vector{Datum}, val::Float64)

    Vector{Datum}(
        [
            d1 + val
            for d1 in d1v
        ]
    )
end

function /(d1v::Vector{Datum}, val::Float64)

    Vector{Datum}(
        [
            d1 / val
            for d1 in d1v
        ]
    )
end

function *(d1v::Vector{Datum}, val::Float64)

    Vector{Datum}(
        [
            d1 * val
            for d1 in d1v
        ]
    )
end



# Data
struct Data
    d::Vector{Datum}

    function Data(
        dd::Vector{Datum}
    )
        new(dd)
    end

    function Data(
        arr::Vector{Tuple{Float64, Float64}}
    )

        d_vec = [
            Datum(x, val) for (x, val) in arr
        ]

        new(d_vec)
    end

    function Data(
        radius_data::Vector{Float64},
        vals::Vector{Float64}
    )
        d_vec = [
            Datum(x, val)
            for (x, val) in zip(radius_data, vals)
        ]

        new(d_vec)
    end

    function Data()
        new([])
    end
end

function -(ds::Data)
    Data(
        [
            (d.x, -d.val) for d in ds.d
        ]
    )
end

function -(ds1::Data, ds2::Data)
    Data(
        ds1.d - ds2.d
    )
end

function -(data::Data, val::Float64)
    Data(
        data.d - val
    )
end

function -(data::Data, val::Int64)
    data - float(val)
end

function +(ds1::Data, ds2::Data)
    Data(
        ds1.d + ds2.d
    )
end

function +(data::Data, val::Float64)
    Data(
        data.d + val
    )
end

function +(data::Data, val::Int64)
    data + float(val)
end

function /(data::Data, denom::Float64)
    Data(
        data.d / denom
    )
end

function /(data::Data, denom::Int64)
    data / float(denom)
end

function *(data::Data, denom::Float64)
    Data(
        data.d * denom
    )
end

function *(data::Data, denom::Int64)
    data * float(denom)
end

function *(denom::Float64, data::Data)
    data * denom
end

function *(denom::Int64, data::Data)
    data * float(denom)
end

function push!(ds::Data, d::Datum)
    push!(ds.d, d)
end


