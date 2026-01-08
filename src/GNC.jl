module GNC

using ..BaseDefs
using StaticArrays

struct GyroSensor
    r::SVector{Float64, 3}
    ψ::Float64
    θ::Float64
    ϕ::Float64
    σ::SVector{Float64}
    b::SVector{Float64}
    T::Float64
end

function takemeasure(sv, s::GyroSensor)

end

struct RollControl

end

"""
    estimate(x̂ₖ₋₁, yₖ, uₖ₋₁)

Determine state estimate `x̂ₖ` given current state measure `yₖ` and previous state estimate `xₖ₋₁`
and applied control `uₖ₋₁`.

# Inputs:
- x̂ₖ₋₁: state estimate at instant `k - 1`.
- yₖ: state measure at instant `k`.
- uₖ₋₁: applied control at instant `k - 1`.
"""
function estimate(x̂ₖ₋₁, yₖ, uₖ₋₁)

end

function control()

end

end