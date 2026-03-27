module Navigation

using StaticArrays

export IMUSensor, takemeasure

struct IMUSensor
    r::SVector{3, Float64}
    σ_gyro::Float64
    σ_accl::Float64
    bias_gyro::SVector{3, Float64}
    bias_accl::SVector{3, Float64}
end

function from_dict(dict::AbstractDict)
    return IMUSensor(
        SVector{3}(dict["position"]),
        dict["sigma_gyro"],
        dict["sigma_accl"],
        dict["bias_gyro"],
        dict["bias_accl"],
    )
end

end