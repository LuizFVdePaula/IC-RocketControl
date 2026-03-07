module Simulate

export simulate, postprocess

using ..BaseDefs: calc_ϕ, calc_θ, calc_ψ, rotXYZ
using ..StageDefs: Stage
using ..EnvironmentDefs
using ..Dynamics
using ..GNC
using ..RK4Solver
using ControlSystems
using DataFrames
using LinearAlgebra
using StaticArrays

"""
Given sample time `Ts`
Initialize current state vector `x`
Initialize observer `x̂`
Initialize control `u`
Loop:
- Current time instant `tₖ`
- Take measure `yₖ` based on `xₖ`
- Update `x̂ₖ` based on `x̂ₖ₋₁`, `yₖ` and `uₖ₋₁`
- Calculate control `uₖ` based on `x̂ₖ` and `tₖ`
- Simulate from `tₖ` to `tₖ₊₁ = tₖ + Ts`
- Update `xₖ₊₁` to the last simulated instant of time `tₖ₊₁`
"""
function simulate(stg::Stage, env::Environment, sv₀, trange::AbstractRange, method::KalmanMethod)
    model = DynamicModel(stg)
    imu = IMUSensor(SVector(-0.4, 0, 0), deg2rad(0.2), 0.2, zeros(SVector{3}), zeros(SVector{3}))

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
    u = zeros(3)
    P = method.P₀
    for (i, t) ∈ enumerate(trange[begin:end-1])
        sysc = continuousmodel(model, x̂, sv[8], -sv[3], t)
        sysd = c2d(sysc, Ts)

        y = takemeasure(imu, sv, u, stg, env, t)
        x̂, P = estimate(x̂, y, u, sysd, P, method)
        u = control(x̂, sysd)
        sol = solve(sv, u, stg, env, range(t, t + Ts, step = dt))

        sv_historic[:, 1+(i-1)*n:1+i*n] .= sol
        u_historic[:, 1+(i-1)*n:1+i*n] .= u
        x̂_historic[:, 1+(i-1)*n:1+i*n] .= x̂
        y_historic[:, 1+(i-1)*n:1+i*n] .= y
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
    δp = sv_historic[14, :]
    δq = sv_historic[15, :]
    δr = sv_historic[16, :]
    p̂ = x̂_historic[1, :]
    q̂ = x̂_historic[2, :]
    r̂ = x̂_historic[3, :]
    α̂ = x̂_historic[4, :]
    β̂ = x̂_historic[5, :]
    p_meas = y_historic[1, :]
    q_meas = y_historic[2, :]
    r_meas = y_historic[3, :]
    ẇ_meas = y_historic[4, :]
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
    up = u_historic[1, :]
    uq = u_historic[2, :]
    ur = u_historic[3, :]
    u̇ = [diff(u); 0] / step(ts)
    v̇ = [diff(v); 0] / step(ts)
    ẇ = [diff(w); 0] / step(ts)
    df = DataFrame(
        "x" => x, "y" => y, "h" => h,
        "q0" => q0, "q1" => q1, "q2" => q2, "q3" => q3,
        "u" => u, "v" => v, "w" => w, "du" => u̇, "dv" => v̇, "dw" => ẇ,
        "p" => p, "q" => q, "r" => r,
        "ϕ" => ϕ, "θ" => θ, "ψ" => ψ,
        "αT" => αT, "ϕA" => ϕA, "α" => α, "β" => β, "rwla" => rwla,
        "up" => up, "uq" => uq, "ur" => ur, "δp" => δp, "δq" => δq, "δr" => δr,
        "p̂" => p̂, "q̂" => q̂, "α̂" => α̂, "r̂" => r̂, "β̂" => β̂,
        "p_meas" => p_meas, "q_meas" => q_meas, "dw_meas" => ẇ_meas, "r_meas" => r_meas, "dv_meas" => v̇_meas
    )
    return df
end

end