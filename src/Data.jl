

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

function -(ds::Data)
    Data(
        [
            (d.x, -d.val) for d in ds.d
        ]
    )
end
