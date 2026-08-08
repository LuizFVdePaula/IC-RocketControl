module GNC

using ..BaseDefs
using ..StageDefs: Stage
using ..EnvironmentDefs
using ..Aerodynamics
using ..Dynamics
using ..Navigation: ESKFState, predict, update_baro, update_pos, update_vel, update_rate
using ControlSystems
using ISAtmosphere
using LinearAlgebra
using StaticArrays

export DynamicModel, ControlParameters, takemeasure, estimate, control

struct DynamicModel
    Lref::Float64
    Sref::Float64
    XR::Float64
    τ::Float64
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
    return DynamicModel(
        stg.aed.Lref,
        stg.aed.Sref,
        stg.aed.XR,
        stg.aed.τ,
        t -> stage_mass(stg, t),
        t -> calc_xcm(stg, t),
        t -> getindex(calc_J(stg, t), 1, 1),
        t -> getindex(calc_J(stg, t), 2, 2),
        stg.aed.base.Clp(M, 0),
        stg.aed.deflection.Clδp(M, 0),
        stg.aed.base.Cmq(M, 0),
        (stg.aed.base.Cm0(M, dα) - stg.aed.base.Cm0(M, 0)) / dα,
        (stg.aed.base.CN0(M, dα) - stg.aed.base.CN0(M, 0)) / dα,
        stg.aed.deflection.Cmδq(M, 0),
        stg.aed.deflection.CNδq(M, 0)
    )
end

struct ControlParameters{T}
    Ts_imu::T
    Ts_baro::T
    # ESKF Parameters
    Q_eskf::SMatrix{15, 15, T, 225}
    R_baro::T
    
    # Attitude Autopilot Parameters
    ωn_pitch::T
    ζ_pitch::T
    ωn_roll::T
    ζ_roll::T
    
    # Energy Management Target
    target_apogee::T
end

function takemeasure(sv, u, stg::Stage, env::Environment, t)
    dsv = dynamics(sv, u, stg, env, t)
    TBG = rotXYZ(sv[4], sv[5], sv[6], sv[7])
    g = TBG * SVector(0, 0, gravity)
    ω = SVector{3}(sv[11:13])
    ω̇ = SVector{3}(dsv[11:13])
    acm = SVector{3}(dsv[8:10])
    xcm = calc_xcm(stg, t)
    ρ⃗ = stg.imu.r - SVector(xcm, 0, 0)
    
    # True acceleration at IMU location
    a_true = acm - g + ω̇ × ρ⃗ + ω × (ω × ρ⃗)
    
    # Sensor Noise
    imu_accel = a_true + stg.imu.σ_accl * randn(SVector{3, Float64}) + stg.imu.bias_accl
    imu_gyro = ω + stg.imu.σ_gyro * randn(SVector{3, Float64}) + stg.imu.bias_gyro
    
    # Barometer (simulate altitude measurement)
    h_true = -sv[3]
    baro_h = h_true + 1.0 * randn() # 1m std dev
    
    return imu_accel, imu_gyro, baro_h
end

function estimate(eskf_state::ESKFState, imu_accel, imu_gyro, baro_h, t, is_baro_tick, cp::ControlParameters)
    # Predict using non-linear kinematics at IMU rate
    state_pred = predict(eskf_state, imu_accel, imu_gyro, cp.Ts_imu, cp.Q_eskf)
    
    if t < 0.0
        # Pre-launch calibration on the rail
        # Very high confidence measurements (low noise) for zero states
        R_pos = SMatrix{3, 3, Float64}(I(3) * 1e-6)
        R_vel = SMatrix{3, 3, Float64}(I(3) * 1e-6)
        R_rate = SMatrix{3, 3, Float64}(I(3) * 1e-6)
        
        # We know position is [0, 0, -1.5], velocity is 0, angular rate is 0
        p_meas = SVector(0.0, 0.0, -1.5)
        v_meas = SVector(0.0, 0.0, 0.0)
        w_meas = SVector(0.0, 0.0, 0.0)
        
        state_upd = update_pos(state_pred, p_meas, R_pos)
        state_upd = update_vel(state_upd, v_meas, R_vel)
        state_upd = update_rate(state_upd, w_meas, imu_gyro, R_rate)
        return state_upd
    else
        if is_baro_tick
            # In-flight update using Barometer altitude
            state_upd = update_baro(state_pred, baro_h, cp.R_baro)
            return state_upd
        else
            return state_pred
        end
    end
end

function control(eskf_state::ESKFState, imu_gyro, dm::DynamicModel, t, cp::ControlParameters)
    # 1. Extract Navigation Data
    # Angular rates (corrected for estimated bias)
    p_est, q_est, r_est = imu_gyro - eskf_state.wb
    
    # Attitude Angles (extract from ESKF quaternion)
    ϕ = calc_ϕ(eskf_state.q...)
    θ = calc_θ(eskf_state.q...)
    ψ = calc_ψ(eskf_state.q...)
    
    # 2. Attitude Error (Target is vertical ascent: θ=80°, ϕ=0, ψ=0)
    err_ϕ = 0.0 - ϕ
    err_θ = deg2rad(80) - θ
    err_ψ = 0.0 - ψ
    
    # 3. Dynamic scheduling (Calculate moments effectiveness online)
    h_est = -eskf_state.p[3]
    V_est = max(norm(eskf_state.v), 1.0) # Clamp to avoid div by zero
    ρ = calc_ρ(h_est)
    q_bar = 0.5 * ρ * V_est^2
    
    # Pitch Control Effectiveness
    M_delta_q = q_bar * dm.Sref * dm.Lref * dm.Cmδq / dm.Jyy(t)
    M_q = q_bar * dm.Sref * dm.Lref^2 * dm.Cmq / (2 * V_est * dm.Jyy(t))
    
    # Roll Control Effectiveness
    M_delta_p = q_bar * dm.Sref * dm.Lref * dm.Clδp / dm.Jxx(t)
    M_p = q_bar * dm.Sref * dm.Lref^2 * dm.Clp / (2 * V_est * dm.Jxx(t))
    
    # 4. Successive Loop Closure / PD (Adaptive Gains)
    sign_Mdq = sign(M_delta_q) == 0 ? -1.0 : sign(M_delta_q)
    Mdq_safe = sign_Mdq * max(abs(M_delta_q), 1e-3)
    
    sign_Mdp = sign(M_delta_p) == 0 ? -1.0 : sign(M_delta_p)
    Mdp_safe = sign_Mdp * max(abs(M_delta_p), 1e-3)

    # Calculate required gains for target natural frequency and damping
    Kp_pitch = cp.ωn_pitch^2 / Mdq_safe
    Kd_pitch = (2 * cp.ζ_pitch * cp.ωn_pitch + M_q) / Mdq_safe
    
    Kp_roll = cp.ωn_roll^2 / Mdp_safe
    Kd_roll = (2 * cp.ζ_roll * cp.ωn_roll + M_p) / Mdp_safe
    
    # PD Control laws
    δq_cmd = Kp_pitch * err_θ - Kd_pitch * q_est
    δr_cmd = Kp_pitch * err_ψ - Kd_pitch * r_est # Yaw uses pitch dynamics (symmetric)
    δp_cmd = Kp_roll * err_ϕ - Kd_roll * p_est
    
    # 5. Apogee Control (Energy Management via Symmetric Drag)
    δ_d = 0.0
    if t > 3.0 # Wait for engine burnout (e.g. tb = 3s)
        # Kinematic apogee prediction
        h_apogee = h_est + eskf_state.v[3]^2 / (2 * gravity) 
        
        # If tracking too high, symmetrically deflect fins to increase drag!
        if h_apogee > cp.target_apogee && eskf_state.v[3] < 0 # Moving up
            δ_d = deg2rad(5) # Calculate 5 degrees of drag
        end
    end
    
    # 6. Final Deflections
    # We don't apply control on the rail (t < 1.0 or h < 5.0)
    f_rail = t < 1.5 ? 0.0 : 1.0
    
    u_p = clamp(δp_cmd, deg2rad(-10), deg2rad(10)) * f_rail
    
    # Note: δ_d is calculated but commented out until the 4-fin bijection mapping 
    # and aerodynamic coefficients are updated to accept a pure drag command.
    u_q = clamp(δq_cmd, deg2rad(-10), deg2rad(10)) * f_rail # + δ_d
    u_r = clamp(δr_cmd, deg2rad(-10), deg2rad(10)) * f_rail # + δ_d
    
    return [u_p, u_q, u_r]
end

end