module Simulate

export simulate, postprocess

using ..BaseDefs
using ..StageDefs: Stage
using ..EnvironmentDefs
using ..Dynamics
using ..Navigation: ESKFState
using ..GNC
using ..RK4Solver
using DataFrames
using LinearAlgebra
using StaticArrays

"""
    simulate(stg::Stage, env::Environment, cp::ControlParameters)

Run a full simulation integrating the true dynamics (RK4) while closing the loop 
at period `Ts` with an ESKF estimator and an Attitude Autopilot.
"""
function simulate(stg::Stage, env::Environment, cp::ControlParameters; ctrl_func=control_fpa)
    # Initial conditions
    ϕ0 = deg2rad(0)
    θ0 = deg2rad(80)
    ψ0 = 0
    q0₀ = calc_q0(ϕ0, θ0, ψ0)
    q1₀ = calc_q1(ϕ0, θ0, ψ0)
    q2₀ = calc_q2(ϕ0, θ0, ψ0)
    q3₀ = calc_q3(ϕ0, θ0, ψ0)
    
    # State Vector: [x, y, z, q0, q1, q2, q3, u, v, w, p, q, r, δp, δq, δr, δpdot, δqdot, δrdot]
    sv₀ = [0.0, 0.0, 0.0, q0₀, q1₀, q2₀, q3₀, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    model = DynamicModel(stg)

    # Simulation step sizes
    n = 1 # Run exactly 1 physical timestep between IMU samples
    Ts_imu = cp.Ts_imu
    trange = range(-10.0, 20.0, step = Ts_imu)
    dt = Ts_imu / n
    N = length(range(trange[begin], trange[end]; step = dt))
    
    # Historic arrays
    sv_historic = zeros(length(sv₀), N)
    x̂_historic = zeros(16, N) # [p(3), v(3), q(4), ab(3), wb(3)]
    u_historic = zeros(3, N)  # [u_p, u_q, u_r]
    y_historic = zeros(7, N)  # [accel(3), gyro(3), baro(1)]
    
    sv = sv₀
    u = zeros(3)
    
    # Initialize ESKF
    P0 = SMatrix{15, 15, Float64}(diagm([
        fill(0.1^2, 3);           # Position
        fill(0.01^2, 3);          # Velocity
        fill(deg2rad(3)^2, 3);    # Attitude
        fill(0.1^2, 3);          # Accel bias
        fill(deg2rad(0.1)^2, 3)   # Gyro bias
    ]))

    eskf_state = ESKFState(
        SVector{3}(0.0, 0.0, 0.0),
        SVector{3}(0.0, 0.0, 0.0),
        SVector{4}(q0₀, q1₀, q2₀, q3₀),
        SVector{3}(0.0, 0.0, 0.0),
        SVector{3}(0.0, 0.0, 0.0),
        P0
    )
    
    ticks_per_baro = round(Int, cp.Ts_baro / cp.Ts_imu)
    ticks_per_control = round(Int, cp.Ts_control / cp.Ts_imu)
    baro_h_hold = -sv[3]

    for (i, t) ∈ enumerate(trange[begin:end-1])
        # 1. Take Measurement (IMU and Barometer)
        imu_accel, imu_gyro, baro_h_new = takemeasure(sv, u, stg, env, t)
        
        # Check if it's time for a barometer tick
        is_baro_tick = (i - 1) % ticks_per_baro == 0
        if is_baro_tick
            baro_h_hold = baro_h_new
        end
        
        # 2. ESKF Estimate
        eskf_state = estimate(eskf_state, imu_accel, imu_gyro, baro_h_hold, t, is_baro_tick, cp)
        
        # 3. Control (runs at Control rate)
        is_control_tick = (i - 1) % ticks_per_control == 0
        if is_control_tick
            if isnothing(ctrl_func)
                u = [0.0, 0.0, 0.0]
            else
                u = ctrl_func(eskf_state, imu_gyro, model, t, cp)
            end
        end
        
        # 4. Simulate physics for this step
        sv = solve!(view(sv_historic, :, 1+(i-1)*n:1+i*n), sv, u, stg, env, range(t, t + Ts_imu, step = dt))

        # 5. Log variables
        u_historic[:, 1+(i-1)*n:1+i*n] .= u
        x̂_historic[1:3, 1+(i-1)*n:1+i*n] .= eskf_state.p
        x̂_historic[4:6, 1+(i-1)*n:1+i*n] .= eskf_state.v
        x̂_historic[7:10, 1+(i-1)*n:1+i*n] .= eskf_state.q
        x̂_historic[11:13, 1+(i-1)*n:1+i*n] .= eskf_state.ab
        x̂_historic[14:16, 1+(i-1)*n:1+i*n] .= eskf_state.wb
        
        y_historic[1:3, 1+(i-1)*n:1+i*n] .= imu_accel
        y_historic[4:6, 1+(i-1)*n:1+i*n] .= imu_gyro
        y_historic[7, 1+(i-1)*n:1+i*n] .= baro_h_hold
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
    
    # Observer states (ESKF)
    p̂_ned = x̂_historic[1:3, :]
    v̂_ned = x̂_historic[4:6, :]
    q̂_0 = x̂_historic[7, :]
    q̂_1 = x̂_historic[8, :]
    q̂_2 = x̂_historic[9, :]
    q̂_3 = x̂_historic[10, :]
    ϕ̂ = @. calc_ϕ(q̂_0, q̂_1, q̂_2, q̂_3)
    θ̂ = @. calc_θ(q̂_0, q̂_1, q̂_2, q̂_3)
    ψ̂ = @. calc_ψ(q̂_0, q̂_1, q̂_2, q̂_3)
    a_bias = x̂_historic[11:13, :]
    w_bias = x̂_historic[14:16, :]
    
    # Measurements
    imu_accel_x = y_historic[1, :]
    imu_accel_y = y_historic[2, :]
    imu_accel_z = y_historic[3, :]
    imu_gyro_p = y_historic[4, :]
    imu_gyro_q = y_historic[5, :]
    imu_gyro_r = y_historic[6, :]
    baro_h = y_historic[7, :]
    
    # True aerodyanmics
    αT = zeros(length(ts))
    ϕA = zeros(length(ts))
    γ = zeros(length(ts))
    vN = zeros(length(ts))
    vE = zeros(length(ts))
    vD = zeros(length(ts))
    for (i, sv) in enumerate(eachcol(sv_historic))
        local h = -sv[3]
        local V = sv[8:10]
        local TBG = rotXYZ(sv[4], sv[5], sv[6], sv[7])
        local vBA = V - TBG * windspeed(env, h)
        local v_NED = TBG' * V
        αT[i] = calc_αT(vBA)
        ϕA[i] = calc_ϕA(vBA)
        γ[i] = atan(-v_NED[3], sqrt(v_NED[1]^2 + v_NED[2]^2))
        vN[i] = v_NED[1]
        vE[i] = v_NED[2]
        vD[i] = v_NED[3]
    end
    α = @. atan(tan(αT) * cos(ϕA))
    β = @. asin(sin(αT) * sin(ϕA))
    
    up = u_historic[1, :]
    uq = u_historic[2, :]
    ur = u_historic[3, :]
    
    df = DataFrame(
        "x" => x, "y" => y, "h" => h, "q0" => q0, "q1" => q1, "q2" => q2, "q3" => q3,
        "u" => u, "v" => v, "w" => w, "p" => p, "q" => q, "r" => r, "ϕ" => ϕ, "θ" => θ, "ψ" => ψ,
        "vN" => vN, "vE" => vE, "vD" => vD, "αT" => αT, "ϕA" => ϕA, "α" => α, "β" => β, "gamma" => γ,
        "xobs" => p̂_ned[1, :], "yobs" => p̂_ned[2, :], "hobs" => -p̂_ned[3, :],
        "ϕobs" => ϕ̂, "θobs" => θ̂, "ψobs" => ψ̂,
        "vNobs" => v̂_ned[1, :], "vEobs" => v̂_ned[2, :], "vDobs" => v̂_ned[3, :],
        "abx" => a_bias[1, :], "aby" => a_bias[2, :], "abz" => a_bias[3, :],
        "wbx" => w_bias[1, :], "wby" => w_bias[2, :], "wbz" => w_bias[3, :],
        "up" => up, "uq" => uq, "ur" => ur, "δp" => δp, "δq" => δq, "δr" => δr,
        "imu_ax" => imu_accel_x, "imu_ay" => imu_accel_y, "imu_az" => imu_accel_z,
        "imu_p" => imu_gyro_p, "imu_q" => imu_gyro_q, "imu_r" => imu_gyro_r, "baro_h" => baro_h
    )
    return df
end

end