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
    Ac, Bc, Cc, Dc = model_ABCD(h, V, stg, t)
    sysc = ss(Ac, Bc, Cc, Dc)
    sysd = c2d(sysc, 0.1)
    Q = diagm([1 / deg2rad(1)^2, 1 / deg2rad(1)^2, 1 / deg2rad(1)^2, 1 / deg2rad(1)^2, 1 / deg2rad(1)^2])
    R = diagm([1 / deg2rad(10)^2, 1 / deg2rad(10)^2, 1 / deg2rad(10)^2])
    K = lqr(sysd, Q, R)
    u = -K * x̂
    δp = clamp(u[1], deg2rad(-15), deg2rad(15))
    δq = clamp(u[2], deg2rad(-15), deg2rad(15))
    δr = clamp(u[3], deg2rad(-15), deg2rad(15))
    return [δp, δq, δr]
end

function takemeasure(sv, u, stg, env, t)
    # y = [p, q, ẇ, r, v̇]
    dsv = dynamics(sv, u, stg, env, t)
    TBG = rotXYZ(sv[4], sv[5], sv[6], sv[7])
    g = TBG * SVector(0, 0, gravity)
    ω = sv[11:13]
    ω̇ = dsv[11:13]
    acm = dsv[8:10]
    xcm = calc_xcm(stg, t)
    ρ⃗ = SVector(-0.4, 0, 0) - SVector(xcm, 0, 0)
    as = acm - g + ω̇ × ρ⃗ + ω × (ω × ρ⃗) + deg2rad(0.1) * SVector{3}(randn(3, 1))
    ωs = ω + deg2rad(0.1) * SVector{3}(randn(3, 1))
    #p = sv[11] + deg2rad(0.1) * randn()
    #q = sv[12] + deg2rad(0.1) * randn()
    #r = sv[13] + deg2rad(0.1) * randn()
    #v̇ = dsv[9] + 1e-1 * randn() - g[2]
    #ẇ = dsv[10] + 1e-1 * randn() - g[3]
    #y = [p, q, ẇ, r, v̇]
    y = [ωs[1], ωs[2], as[3], ωs[3], as[2]]
    return y
end

"""
    estimate(x̂ₖ, yₖ, uₖ, h, V, stg::Stage, t)

Estimate state `x̂ₖ₊₁` given previous state estimate `x̂ₖ`, state measure `yₖ` and applied control `uₖ`.

# Inputs:
- x̂ₖ: state estimate at instant `k`.
- yₖ: state measure at instant `k`.
- uₖ: applied control at instant `k`.
"""
function estimate(x̂ₖ, yₖ, uₖ, h, V, stg::Stage, t)
    # x = [p, q, α, r, β]
    Ac, Bc, H, D = model_ABCD(h, V, stg, t)
    sysc = ss(Ac, Bc, H, D)
    sysd = c2d(sysc, 0.1)
    A = sysd.A
    B = sysd.B
    σw = 2 / V
    Gα = Ac[2:3, 3] # [Mα; Nα / V]
    Gβ = Ac[4:5, 5] # [Mβ; Nβ / V]
    G = zeros(5, 3)
    G[1, 1] = 5.6
    G[2:3, 2] .= Gα
    G[4:5, 3] .= Gβ
    Σw = diagm([0.02^2, 4σw^2, σw^2])
    Qc = G * Σw * transpose(G)
    R1 = c2d(sysd, Qc, opt = :o)
    R2 = diagm([deg2rad(0.1)^2, deg2rad(0.1)^2, 1e-2, deg2rad(0.1)^2, 1e-2])
    #K = kalman(sysd, R1, R2)
    K = place(A, H, 1e-1 * [1, 2, 1im, -2, -1], :o)
    x̂ₖ₊₁ = (A - K * H) * x̂ₖ + (B - K * D) * uₖ + K * yₖ
    return x̂ₖ₊₁
end

function model_ABCD(h, V, stg::Stage, t)
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
    C = [1 0 0 0 0; 0 1 0 0 0; 0 V Nα 0 0; 0 0 0 1 0; 0 0 0 -V Nβ]
    D = [0 0 0; 0 0 0; 0 Nδq 0; 0 0 0; 0 0 Nδr]
    return (A, B, C, D)
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
- Take measure `yₖ` based on `xₖ`
- Simulate from `tₖ` to `tₖ₊₁ = tₖ + Ts`
- Update `x̂ₖ₊₁` based on `x̂ₖ`, `yₖ` and `uₖ`
- Update `xₖ₊₁` to the last simulated instant of time `tₖ₊₁`
"""
function simulate(stg::Stage, env::Environment, sv₀, trange::AbstractRange)
    model = DynamicModel(stg)
    imu = IMUSensor(SVector(-0.4, 0, 0), deg2rad(1), deg2rad(1), zeros(SVector, 3), zeros(SVector, 3))

    n = 10
    Ts = step(trange)
    dt = Ts / n
    N = length(range(trange[begin], trange[end]; step = dt))
    sv_historic = zeros(length(sv₀), N)
    x̂_historic = zeros(5, N)
    u_historic = zeros(3, N)
    y_historic = zeros(5, N)
    sv = sv₀
    x̂ = zeros(5)
    u = calc_control(x̂, t, -sv[3], sv[8], stg)
    for (i, t) ∈ enumerate(trange[begin:end-1])
        sysc = continuousmodel(model, x̂, sv[8], -sv[3])
        sysd = c2d(sysc, Ts)

        y = takemeasure(imu, sv, u, stg, env, t)
        u = control(x̂, sysd)
        sol = solve(sv, u, stg, env, range(t, t + Ts, step = dt))
        sv_historic[:, 1+(i-1)*n:1+i*n] .= sol
        u_historic[:, 1+(i-1)*n:1+i*n] .= u
        x̂_historic[:, 1+(i-1)*n:1+i*n] .= x̂
        y_historic[:, 1+(i-1)*n:1+i*n] .= y
        x̂ = estimate(x̂, y, u, sysc, sysd)
        sv = sol[:, end]
    end
    return sv_historic, u_historic, x̂_historic, y_historic
end

function postprocess(stg, env, sv_historic, u_historic, x̂_historic, y_historic, ts)
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
    p_meas = y_historic[1, :]
    q_meas = y_historic[2, :]
    ẇ_meas = y_historic[3, :]
    r_meas = y_historic[4, :]
    v̇_meas = y_historic[5, :]
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
    u̇ = [diff(u); 0] / step(ts)
    v̇ = [diff(v); 0] / step(ts)
    ẇ = [diff(w); 0] / step(ts)
    df = DataFrame(
        "x" => x, "y" => y, "h" => h,
        "q0" => q0, "q1" => q1, "q2" => q2, "q3" => q3,
        "u" => u, "v" => v, "w" => w, "du" => u̇, "dv" => v̇, "dw" => ẇ,
        "p" => p, "q" => q, "r" => r,
        "ϕ" => ϕ, "θ" => θ, "ψ" => ψ,
        "αT" => αT, "ϕA" => ϕA, "α" => α, "β" => β,
        "rwla" => rwla,
        "δp" => δp, "δq" => δq, "δr" => δr,
        "p̂" => p̂, "q̂" => q̂, "α̂" => α̂, "r̂" => r̂, "β̂" => β̂,
        "p_meas" => p_meas, "q_meas" => q_meas, "dw_meas" => ẇ_meas, "r_meas" => r_meas, "dv_meas" => v̇_meas
    )
    return df
end

end