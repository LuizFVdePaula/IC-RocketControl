module Simulate

export simulate, postprocess

using ..BaseDefs: calc_ϕ, calc_θ, calc_ψ, rotXYZ
using ..StageDefs: Stage
using ..Aerodynamics
using ..EnvironmentDefs
using ..Dynamics
using ..RK4Solver
using ISAtmosphere
using ControlSystems
using DataFrames
using LinearAlgebra
using StaticArrays

function calc_control(x̂, t, h, V, stg)
    # x̂ = [p, q, α, r, β]
    Ac, Bc = model_AB(h, V, stg, t)
    H = [1 0 0 0 0; 0 1 0 0 0; 0 0 0 1 0]
    sysc = ss(Ac, Bc, H, 0)
    sysd = c2d(sysc, 0.1)
    Q = diagm([1 / deg2rad(1)^2, 1 / deg2rad(2)^2, 1 / deg2rad(1)^2, 1 / deg2rad(1)^2, 1 / deg2rad(1)^2])
    R = diagm([1 / deg2rad(5)^2, 1 / deg2rad(5)^2, 1 / deg2rad(5)^2])
    K = lqr(sysd, Q, R)
    u = -K * x̂
    δp = clamp(u[1], -deg2rad(15), deg2rad(15))
    δq = clamp(u[2], -deg2rad(15), deg2rad(15))
    δr = clamp(u[3], -deg2rad(15), deg2rad(15))
    return [δp, δq, δr]
end

function takemeasure(sv)
    # y = [p, q, r]
    p = sv[11] + 1e-3 * randn()
    q = sv[12] + 1e-3 * randn()
    r = sv[13] + 1e-3 * randn()
    y = [p, q, r]
    return y
end

"""
    estimate(x̂ₖ₋₁, yₖ, uₖ₋₁, h, V, stg::Stage, t)

Determine state estimate `x̂ₖ` given current state measure `yₖ` and previous state estimate `xₖ₋₁`
and applied control `uₖ₋₁`.

# Inputs:
- x̂ₖ₋₁: state estimate at instant `k - 1`.
- yₖ: state measure at instant `k`.
- uₖ₋₁: applied control at instant `k - 1`.
"""
function estimate(x̂ₖ₋₁, yₖ, uₖ₋₁, h, V, stg::Stage, t)
    Ac, Bc = model_AB(h, V, stg, t)
    H = [1 0 0 0 0; 0 1 0 0 0; 0 0 0 1 0]
    sysc = ss(Ac, Bc, H, 0)
    sysd = c2d(sysc, 0.1)
    A = sysd.A
    B = sysd.B
    x̂ₖ_pri = A * x̂ₖ₋₁ + B * uₖ₋₁
    ŷₖ_pri = H * x̂ₖ_pri
    R1 = diagm([0.005, 0.03, 0.02, 0.03, 0.02])
    R2 = diagm([1e-3, 1e-3, 1e-3])
    #K = place(A, H, 1e-1 * [1, 2, 3, 4, 5], :o; direct = true)
    K = kalman(sysd, R1, R2; direct = true)
    x̂ₖ_pos = x̂ₖ_pri + K * (yₖ - ŷₖ_pri)
    return x̂ₖ_pos
end

function model_AB(h, V, stg::Stage, t)
    ρ = ρ_kg_m³(p_Pa(h), T_K(h))
    S = stg.aed.Sref
    c = stg.aed.Lref
    m = stage_mass(stg, t)
    J = calc_J(stg, t)
    Ixx = J[1, 1]
    Iyy = J[2, 2]
    k = ρ * V^2 * S / 2
    M = V / a_m_s(T_K(h))
    XCG = calc_xcm(stg, t) / c
    ΔXCG = stg.aed.XR - XCG

    Clp = (getCl(stg.aed, M, 0, 0, 1e-2, 0) - getCl(stg.aed, M, 0, 0, 0, 0)) / 1e-2
    Clδp = (getCl(stg.aed, M, 0, 0, 0, 1e-2) - getCl(stg.aed, M, 0, 0, 0, 0)) / 1e-2
    Cmq = (getCm(stg.aed, M, 0, 0, ΔXCG, 1e-2, 0) - getCm(stg.aed, M, 0, 0, ΔXCG, 0, 0)) / 1e-2
    Cmα = (getCm(stg.aed, M, 1e-2, 0, ΔXCG, 0, 0) - getCm(stg.aed, M, 0, 0, ΔXCG, 0, 0)) / 1e-2
    CNα = (getCN(stg.aed, M, 1e-2, 0, 0) - getCN(stg.aed, M, 0, 0, 0)) / 1e-2
    Cmδq = (getCm(stg.aed, M, 0, 0, ΔXCG, 0, 1e-2) - getCm(stg.aed, M, 0, 0, ΔXCG, 0, 0)) / 1e-2
    CNδq = (getCN(stg.aed, M, 0, 0, 1e-2) - getCN(stg.aed, M, 0, 0, 0)) / 1e-2
    Cnr = Cmq
    Cnβ = -Cmα
    CYβ = CNα
    Cnδr = Cmδq
    CYδr = -CNδq

    Mp = k * c^2 / (2 * V * Ixx) * Clp
    Mδp = k * c / Ixx * Clδp
    Mq = k * c^2 / (2 * V * Iyy) * Cmq
    Mα = k * c / Iyy * Cmα
    Nα = k * CNα / m
    Mδq = k * c / Iyy * Cmδq
    Nδq = k * CNδq / m
    Mr = k * c^2 / (2 * V * Iyy) * Cnr
    Mβ = k * c / Iyy * Cnβ
    Nβ = k * CYβ / m
    Mδr = k * c / Iyy * Cnδr
    Nδr = k * CYδr / m

    Arol = Mp
    Brol = Mδp
    Alon = [Mq Mα; 1 Nα / V]
    Blon = [Mδq; Nδq / V]
    Alat = [Mr Mβ; -1 Nβ / V]
    Blat = [Mδr; Nδr / V]

    A = [Arol zeros(1, 4); zeros(2, 1) Alon zeros(2, 2); zeros(2, 3) Alat]
    B = [Brol zeros(1, 2); zeros(2, 1) Blon zeros(2, 1); zeros(2, 2) Blat]
    return (A, B)
end

"""
1 -> sv₀ (t = 0.0)
solve resolve de [1:11] (t = 0.0 até 1.0)
append sol[2:11] (t = 0.1 até 1.0)
sv = sol[end]

Given sample time `Ts`
Initialize current state vector `x`
Initialize observer `x̂`
Loop:
- Current time instant `tₖ`
- Calculate control `uₖ` based on `x̂ₖ` and `tₖ`
- Simulate from `tₖ` to `tₖ₊₁ = tₖ + Ts`
- Update `xₖ₊₁` to the last simulated instant of time `tₖ₊₁`
- Take measure `yₖ₊₁` based on `xₖ₊₁`
- Update `x̂ₖ₊₁` based on `x̂ₖ`, `yₖ₊₁` and `uₖ`
"""
function simulate(stg::Stage, env::Environment, sv₀, trange::AbstractRange)
    n = 10
    Ts = step(trange)
    dt = Ts / n
    N = length(range(trange[begin], trange[end]; step = dt))
    sv_historic = zeros(length(sv₀), N)
    x̂_historic = zeros(5, N)
    u_historic = zeros(3, N)
    x̂ = zeros(5)
    sv = sv₀
    for (i, t) ∈ enumerate(trange[begin:end-1])
        u = calc_control(x̂, t, -sv[3], sv[8], stg)
        sol = solve(sv, u, stg, env, range(t, t + Ts, step = dt))
        sv_historic[:, 1+(i-1)*n:1+i*n] .= sol
        u_historic[:, 1+(i-1)*n:1+i*n] .= u
        x̂_historic[:, 1+(i-1)*n:1+i*n] .= x̂
        sv = sol[:, end]
        y = takemeasure(sv)
        x̂ = estimate(x̂, y, u, -sv[3], sv[8], stg, t + Ts)
    end
    return sv_historic, u_historic, x̂_historic
end

function postprocess(stg, env, sv_historic, u_historic, x̂_historic)
    x = sv_historic[1, :]
    y = sv_historic[2, :]
    h = -sv_historic[3, :]
    q0 = sv_historic[4, :]
    q1 = sv_historic[5, :]
    q2 = sv_historic[6, :]
    q3 = sv_historic[7, :]
    ϕ = @. calc_ϕ(q0, q1, q2, q3)
    θ = @. calc_θ(q0, q1, q2, q3)
    ψ = @. calc_ψ(q0, q1, q2, q3)
    u = sv_historic[8, :]
    v = sv_historic[9, :]
    w = sv_historic[10, :]
    p = sv_historic[11, :]
    q = sv_historic[12, :]
    r = sv_historic[13, :]
    p̂ = x̂_historic[1, :]
    q̂ = x̂_historic[2, :]
    α̂ = x̂_historic[3, :]
    r̂ = x̂_historic[4, :]
    β̂ = x̂_historic[5, :]
    αT = map(eachcol(sv_historic)) do sv
        local h = -sv[3]
        local V = sv[8:10]
        local TBG = rotXYZ(sv[4], sv[5], sv[6], sv[7])
        local vBA = V - TBG * windspeed(env, h)
        calc_αT(vBA)
    end
    ϕA = map(eachcol(sv_historic)) do sv
        local h = -sv[3]
        local V = sv[8:10]
        local TBG = rotXYZ(sv[4], sv[5], sv[6], sv[7])
        local vBA = V - TBG * windspeed(env, h)
        calc_ϕA(vBA)
    end
    α = @. atan(tan(αT) * cos(ϕA))
    β = @. asin(sin(αT) * sin(ϕA))
    rwla = map(eachcol(sv_historic)) do sv
        local TBG = rotXYZ(sv[4:7]...)
        local yB = transpose(TBG) * [0; 1; 0]
        atan(yB[2], yB[1])
    end
    δp = u_historic[1, :]
    δq = u_historic[2, :]
    δr = u_historic[3, :]
    df = DataFrame(
        "x" => x, "y" => y, "h" => h,
        "q0" => q0, "q1" => q1, "q2" => q2, "q3" => q3,
        "u" => u, "v" => v, "w" => w,
        "p" => p, "q" => q, "r" => r,
        "ϕ" => ϕ, "θ" => θ, "ψ" => ψ,
        "αT" => αT, "ϕA" => ϕA, "α" => α, "β" => β,
        "rwla" => rwla,
        "δp" => δp, "δq" => δq, "δr" => δr,
        "p̂" => p̂, "q̂" => q̂, "α̂" => α̂, "r̂" => r̂, "β̂" => β̂
    )
    return df
end

end