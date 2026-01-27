module GNC

using ..BaseDefs
using ..StageDefs
using ..EnvironmentDefs
using ..Dynamics: dynamics
using LinearAlgebra
using StaticArrays

"""
Complete state: x = [x, y, z, q₀, q₁, q₂, q₃, u, v, w, p, q, r, δp, δq, δr]
Estimated state: x̂ = [p, q, r, α, β]
Measured: z = [p, q, r, v̇, ẇ]
Control input: u = [up, uq, ur]
Control law: u = -K ⋅ x̂
"""

struct IMUSensor
    r::SVector{Float64, 3}
    σ_gyro::Float64
    σ_accl::Float64
    bias_gyro::SVector{Float64, 3}
    bias_accl::SVector{Float64, 3}
end

function takemeasure(imu::IMUSensor, sv, u, stg::Stage, env::Environment, t)
    # z = [p, q, r, v̇, ẇ]
    dsv = dynamics(sv, u, stg, env, t)
    TBG = rotXYZ(sv[4], sv[5], sv[6], sv[7])
    g = TBG * SVector(0, 0, gravity)
    ω = sv[11:13]
    ω̇ = dsv[11:13]
    acm = dsv[8:10]
    xcm = calc_xcm(stg, t)
    ρ⃗ = imu.r - SVector(xcm, 0, 0)
    as = acm - g + ω̇ × ρ⃗ + ω × (ω × ρ⃗) + imu.σ_accl * SVector{3}(randn(3, 1)) + imu.bias_accl
    ωs = ω + imu.σ_gyro * SVector{3}(randn(3, 1)) + imu.bias_gyro
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
    Lref
    Sref
    XR
    m
    xcm
    Jxx
    Jyy
    Clp
    Clδp
    Cmq
    Cmα
    CNα
    Cmδq
    CNδq
end

function DynamicModel(stg::Stage)
    return DynamicModel(
        stg.aed.Lref,
        stg.aed.Sref,
        stg.aed.XR,
        t -> stage_mass(stg, t),
        t -> calc_xcm(stg, t),
        t -> getindex(calc_J(stg, t), 1, 1),
        t -> getindex(calc_J(stg, t), 2, 2),
        Clp = (getCl(stg.aed, M, 0, 0, 1e-2, 0) - getCl(stg.aed, M, 0, 0, 0, 0)) / 1e-2,
        Clδp = (getCl(stg.aed, M, 0, 0, 0, 1e-2) - getCl(stg.aed, M, 0, 0, 0, 0)) / 1e-2,
        Cmq = (getCm(stg.aed, M, 0, 0, ΔXCG, 1e-2, 0) - getCm(stg.aed, M, 0, 0, ΔXCG, 0, 0)) / 1e-2,
        Cmα = (getCm(stg.aed, M, 1e-2, 0, ΔXCG, 0, 0) - getCm(stg.aed, M, 0, 0, ΔXCG, 0, 0)) / 1e-2,
        CNα = (getCN(stg.aed, M, 1e-2, 0, 0) - getCN(stg.aed, M, 0, 0, 0)) / 1e-2,
        Cmδq = (getCm(stg.aed, M, 0, 0, ΔXCG, 0, 1e-2) - getCm(stg.aed, M, 0, 0, ΔXCG, 0, 0)) / 1e-2,
        CNδq = (getCN(stg.aed, M, 0, 0, 1e-2) - getCN(stg.aed, M, 0, 0, 0)) / 1e-2
    )
end

function continuousmodel(dm::DynamicModel, x̂, V, h)
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

    A = [
        Mp 0  0  0      0
        0  Mq 0  Mα     0
        0  1  0  Nα / V 0
        0  0  Mr 0      Mβ
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
        0 Nδr 0
        0 0   Nδq
    ]
    return ss(A, B, C, D)
end

struct KalmanMethod
    P₀
end

"""
    estimate(x̂ₖ, yₖ, uₖ, sysc, model)

Estimate state `x̂ₖ₊₁` given previous state estimate `x̂ₖ`, state measure `yₖ` and applied control `uₖ`.

# Inputs:
- `x̂ₖ`: state estimate at instant `k`.
- `yₖ`: state measure at instant `k`.
- `uₖ`: applied control at instant `k`.
- `sysc`: continuous-time system model.
- `method`: method (Kalman, Luenberger, Reduced).
"""
function estimate(x̂ₖ, yₖ, uₖ, sysc, sysd)
    # x̂ = [p, q, r, α, β]
    method = "kalman"
    A = sysd.A
    B = sysd.B
    H = sysd.C
    D = sysd.D

    if method == "kalman"
        Ac = sysc.A
        σw = 2 / 100
        Mα, Nα_V = Ac[2:3, 4] # [Mα; Nα / V]
        Mβ, Nβ_V = Ac[4:5, 5] # [Mβ; Nβ / V]
        G = [
            5.6 0    0
            0   Mα   0
            0   0    Mβ
            0   Nα_V 0
            0   0    Nβ_V
        ]
        Σw = diagm([0.02^2, 4σw^2, σw^2])
        Qc = G * Σw * transpose(G)
        R1 = c2d(sysd, Qc, opt = :o)
        R2 = diagm([deg2rad(1)^2, deg2rad(1)^2, deg2rad(1)^2, 1, 1])
        K = kalman(sysd, R1, R2)
    elseif method == "pole"
        K = place(A, H, 1e-1 * [1, 2, 1im, -2, -1], :o)
    else
        error("Method $method not implemented.")
    end

    x̂ₖ₊₁ = (A - K * H) * x̂ₖ + (B - K * D) * uₖ + K * yₖ
    return x̂ₖ₊₁
end

function control(x̂, sysd)
    # x̂ = [p, q, r, α, β]
    Q = diagm([1 / deg2rad(1)^2, 1 / deg2rad(1)^2, 1 / deg2rad(1)^2, 1 / deg2rad(1)^2, 1 / deg2rad(1)^2])
    R = diagm([1 / deg2rad(10)^2, 1 / deg2rad(10)^2, 1 / deg2rad(10)^2])
    L = lqr(sysd, Q, R)
    u = -L * x̂
    δp = clamp(u[1], deg2rad(-15), deg2rad(15))
    δq = clamp(u[2], deg2rad(-15), deg2rad(15))
    δr = clamp(u[3], deg2rad(-15), deg2rad(15))
    return [δp, δq, δr]
end

end