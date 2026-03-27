module GNC

using ..BaseDefs
using ..StageDefs: Stage, IMUSensor
using ..EnvironmentDefs
using ..Aerodynamics
using ..Dynamics
using ControlSystems
using ISAtmosphere
using LinearAlgebra
using StaticArrays

export DynamicModel, KalmanMethod, takemeasure, continuousmodel, estimate, control

"""
Complete state: x = [x, y, z, q₀, q₁, q₂, q₃, u, v, w, p, q, r, δp, δq, δr]
Estimated state: x̂ = [p, q, r, α, β]
Measured: z = [p, q, r, v̇, ẇ]
Control input: u = [up, uq, ur]
Control law: u = -K ⋅ x̂
"""

function takemeasure(sv, u, stg::Stage, env::Environment, t)
    # z = [p, q, r, v̇, ẇ]
    dsv = dynamics(sv, u, stg, env, t)
    TBG = rotXYZ(sv[4], sv[5], sv[6], sv[7])
    g = TBG * SVector(0, 0, gravity)
    ω = sv[11:13]
    ω̇ = dsv[11:13]
    acm = dsv[8:10]
    xcm = calc_xcm(stg, t)
    ρ⃗ = stg.imu.r - SVector(xcm, 0, 0)
    as = acm - g + ω̇ × ρ⃗ + ω × (ω × ρ⃗) + stg.imu.σ_accl * SVector{3}(randn(3, 1)) + stg.imu.bias_accl
    ωs = ω + stg.imu.σ_gyro * SVector{3}(randn(3, 1)) + stg.imu.bias_gyro
    z = SVector{5}([ωs; as[2:3]])
    return z
end

"""
    DynamicModel

The model used in order to construct the continuous-time system.

Sref, Lref => constant
m, Jxx, Jyy, xcm => interpolation of time (feedforward)
Coeff => constant?, consider coeff around XR and transfer Cm inside 'continuousmodel'
"""
struct DynamicModel
    Lref::Float64
    Sref::Float64
    XR::Float64
    m
    xcm
    Jxx
    Jyy
    Clp::Float64
    Clδp::Float64
    Cmq::Float64
    Cmα::Float64
    CNα::Float64
    Cmδq::Float64
    CNδq::Float64
end

function DynamicModel(stg::Stage)
    M = 0.25
    dα = deg2rad(3)
    dδ = deg2rad(3)
    dω = deg2rad(2)
    return DynamicModel(
        stg.aed.Lref,
        stg.aed.Sref,
        stg.aed.XR,
        t -> stage_mass(stg, t),
        t -> calc_xcm(stg, t),
        t -> getindex(calc_J(stg, t), 1, 1),
        t -> getindex(calc_J(stg, t), 2, 2),
        (getCl(stg.aed, M, 0, 0, dω, 0) - getCl(stg.aed, M, 0, 0, 0, 0)) / dω,
        (getCl(stg.aed, M, 0, 0, 0, dδ) - getCl(stg.aed, M, 0, 0, 0, 0)) / dδ,
        (getCm(stg.aed, M, 0, 0, 0, dω, 0) - getCm(stg.aed, M, 0, 0, 0, 0, 0)) / dω,
        (getCm(stg.aed, M, dα, 0, 0, 0, 0) - getCm(stg.aed, M, 0, 0, 0, 0, 0)) / dα,
        (getCN(stg.aed, M, dα, 0, 0) - getCN(stg.aed, M, 0, 0, 0)) / dα,
        (getCm(stg.aed, M, 0, 0, 0, 0, dδ) - getCm(stg.aed, M, 0, 0, 0, 0, 0)) / dδ,
        (getCN(stg.aed, M, 0, 0, dδ) - getCN(stg.aed, M, 0, 0, 0)) / dδ
    )
end

function continuousmodel(dm::DynamicModel, x̂, V, h, t)
    ρ = ρ_kg_m³(p_Pa(h), T_K(h))
    S = dm.Sref
    c = dm.Lref
    m = dm.m(t)
    Jxx = dm.Jxx(t)
    Jyy = dm.Jyy(t)
    k = ρ * V^2 * S / 2
    XCG = dm.xcm(t) / c
    ΔXCG = dm.XR - XCG
    Clp = dm.Clp
    Clδp = dm.Clδp
    Cmq = dm.Cmq
    Cmα = dm.Cmα - dm.CNα * ΔXCG
    CNα = dm.CNα
    Cmδq = dm.Cmδq #- dm.CNδq * ΔXCG # review in aerodynamic model
    CNδq = dm.CNδq
    Cnr = Cmq
    Cnβ = -Cmα
    CYβ = CNα
    Cnδr = Cmδq
    CYδr = -CNδq

    Mp = k * c^2 / (2 * V * Jxx) * Clp
    Mδp = k * c / Jxx * Clδp
    Mq = k * c^2 / (2 * V * Jyy) * Cmq
    Mα = k * c / Jyy * Cmα
    Nα = k * CNα / m
    Mδq = k * c / Jyy * Cmδq
    Nδq = k * CNδq / m
    Mr = k * c^2 / (2 * V * Jyy) * Cnr
    Mβ = k * c / Jyy * Cnβ
    Nβ = k * CYβ / m
    Mδr = k * c / Jyy * Cnδr
    Nδr = k * CYδr / m
    #TODO: include actuator 1st order dynamics
    A = [
        Mp 0  0  0      0
        0  Mq 0  Mα     0
        0  0  Mr 0      Mβ
        0  1  0  Nα / V 0
        0  0  -1 0      Nβ / V
    ]
    B = [
        Mδp 0       0
        0   Mδq     0
        0   0       Mδr
        0   Nδq / V 0
        0   0       Nδr / V
    ]
    C = [
        1 0 0  0  0
        0 1 0  0  0
        0 0 1  0  0
        0 0 -V 0  Nβ
        0 V 0  Nα 0
    ]
    D = [
        0 0   0
        0 0   0
        0 0   0
        0 0   Nδr
        0 Nδq 0
    ]
    return ss(A, B, C, D)
end

struct KalmanMethod
    P₀
    Q
    R
end

"""
    estimate(x̂ₖ₋₁, zₖ, uₖ₋₁, sysd, model)

Estimate state `x̂ₖ` given previous state estimate `x̂ₖ₋₁`, applied control `uₖ₋₁` and current state measure `zₖ`.

# Inputs:
- `x̂ₖ₋₁`: previous state estimate at instant `k - 1`.
- `zₖ`: current state measure at instant of estimation `k`.
- `uₖ₋₁`: applied control between instants `k - 1` and `k`.
- `sysd`: discrete-time system model.
"""
function estimate(x̂ₖ₋₁, zₖ, uₖ₋₁, sysd, Pₖ₋₁, method::KalmanMethod)
    # x̂ = [p, q, r, α, β]
    A = sysd.A
    B = sysd.B
    H = sysd.C
    D = sysd.D

    x̂ₖ⁻ = A * x̂ₖ₋₁ + B * uₖ₋₁
    Pₖ⁻ = A * Pₖ₋₁ * A' + method.Q

    ẑₖ⁻ = H * x̂ₖ⁻  + D * uₖ₋₁
    yₖ = zₖ - ẑₖ⁻
    S = H * Pₖ⁻ * H' + method.R
    K = (Pₖ⁻ * H') / S

    #TODO Joseph form

    Pₖ⁺ = (I - K * H) * Pₖ⁻
    x̂ₖ⁺ = x̂ₖ⁻ + K * yₖ
    return x̂ₖ⁺, Pₖ⁺
end

function control(x̂, sysd, t)
    tref = 1.0
    f = min(t / tref, 1.0)
    # x̂ = [p, q, r, α, β]
    Q = diagm([1 / deg2rad(3)^2, 1 / deg2rad(3)^2, 1 / deg2rad(3)^2, 1 / deg2rad(1.0)^2, 1 / deg2rad(1.0)^2])
    R = diagm([1 / deg2rad(30)^2, 1 / deg2rad(30)^2, 1 / deg2rad(30)^2])
    L = lqr(sysd, Q, R)
    u = -L * x̂
    δp = clamp(u[1], deg2rad(-10), deg2rad(10))
    δq = clamp(u[2], deg2rad(-10), deg2rad(10))
    δr = clamp(u[3], deg2rad(-10), deg2rad(10))
    return f * [δp, δq, δr]
end

end