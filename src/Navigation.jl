module Navigation

using StaticArrays
using LinearAlgebra
using ..BaseDefs: rotXYZ
using ..EnvironmentDefs: gravity

export IMUSensor, ESKFState, predict, update_baro, update_rail

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

"""
    ESKFState

Maintains the nominal state and the 15x15 covariance matrix of the error state.
- `p`: Position (NED)
- `v`: Velocity (NED)
- `q`: Attitude Quaternion (Ground to Body, scalar first)
- `ab`: Accelerometer Bias
- `wb`: Gyroscope Bias
- `P`: 15x15 Error State Covariance Matrix
"""
struct ESKFState
    p::SVector{3, Float64}
    v::SVector{3, Float64}
    q::SVector{4, Float64}
    ab::SVector{3, Float64}
    wb::SVector{3, Float64}
    P::SMatrix{15, 15, Float64, 225}
end

function skew(v)
    return SMatrix{3, 3, Float64}([
        0.0  -v[3]  v[2];
        v[3]  0.0  -v[1];
        -v[2] v[1]  0.0
    ])
end

function quat_mult(q1, q2)
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return SVector{4, Float64}(
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
    )
end

"""
    predict(state::ESKFState, imu_accel, imu_gyro, dt, Q)

Propagate the nominal state using non-linear kinematics and propagate the 
error covariance `P` using the linearized error-state Jacobian.
"""
function predict(state::ESKFState, imu_accel, imu_gyro, dt, Q)
    p = state.p
    v = state.v
    q = state.q
    ab = state.ab
    wb = state.wb
    P = state.P

    # True IMU readings (bias removed)
    a_true = imu_accel - ab
    w_true = imu_gyro - wb

    # 1. Nominal Kinematic Propagation
    TBG = rotXYZ(q...) # Ground to Body DCM
    TGB = transpose(TBG) # Body to Ground DCM
    g_vec = SVector(0.0, 0.0, gravity)
    
    a_ground = TGB * a_true + g_vec
    
    p_new = p + v * dt + 0.5 * a_ground * dt^2
    v_new = v + a_ground * dt
    
    wx, wy, wz = w_true
    Ωquat = SMatrix{4, 4, Float64}([
        0.0 -wx -wy -wz;
        wx 0.0 wz -wy;
        wy -wz 0.0 wx;
        wz wy -wx 0.0
    ])
    q_new = q + 0.5 * dt * Ωquat * q
    q_new = q_new / norm(q_new)

    # 2. Error-state Jacobian F (15x15)
    # Order: δp (1:3), δv (4:6), δθ (7:9), δab (10:12), δwb (13:15)
    # Creating F using block matrices for performance
    F_M = @MMatrix zeros(15, 15)
    F_M[1:3, 4:6] = I(3)
    F_M[4:6, 7:9] = -TGB * skew(a_true)
    F_M[4:6, 10:12] = -TGB
    F_M[7:9, 7:9] = -skew(w_true)
    F_M[7:9, 13:15] = -I(3)
    F = SMatrix{15, 15, Float64}(F_M)

    # Discrete transition matrix (Euler integration)
    Phi = I(15) + F * dt
    
    # 3. Covariance Propagation
    P_new = Phi * P * transpose(Phi) + Q * dt

    return ESKFState(p_new, v_new, q_new, ab, wb, P_new)
end

"""
    inject(state::ESKFState, dx::SVector{15})

Injects the estimated error state `dx` into the nominal state to correct it, 
and then zeroes out the error state.
"""
function inject(state::ESKFState, dx::SVector{15, Float64})
    dp = dx[1:3]
    dv = dx[4:6]
    dth = dx[7:9]
    dab = dx[10:12]
    dwb = dx[13:15]

    p_new = state.p + dp
    v_new = state.v + dv
    
    # Attitude error injection (local error)
    # q_new = q_nom ⊗ [1; dth / 2]
    dq = SVector{4, Float64}(1.0, dth[1] / 2, dth[2] / 2, dth[3] / 2)
    q_new = quat_mult(state.q, dq)
    q_new = q_new / norm(q_new)

    ab_new = state.ab + dab
    wb_new = state.wb + dwb

    G = Matrix{Float64}(I, 15, 15)
    G[7:9, 7:9] .= I(3) - 0.5 * skew(dth)
    P_new = G * state.P * transpose(G)

    return ESKFState(p_new, v_new, q_new, ab_new, wb_new, P_new)
end

"""
    update_baro(state::ESKFState, h_meas, R_baro)

Kalman update step using a barometer altitude measurement.
"""
function update_baro(state::ESKFState, h_meas, R_baro)
    # Barometer measures altitude (negative z in NED)
    z_hat = -state.p[3]
    y = h_meas - z_hat

    # Measurement Jacobian H (1x15) mapping error state to altitude error
    # H = SMatrix{1, 15, Float64}([0 0 -1.0 0 0 0 0 0 0 0 0 0 0 0 0])

    # Innovation covariance: H' P H + R
    S = state.P[3, 3] + R_baro
    
    # Kalman Gain: P H' S^-1
    K = -state.P[:, 3] / S

    # Error state calculation
    dx = K * y
    
    # Joseph form covariance update (numerically stable)
    I_KH = MMatrix{15, 15, Float64}(I(15))
    I_KH[:, 3] .+= K
    P_new = I_KH * state.P * transpose(I_KH) + K * R_baro * transpose(K)

    state_updated = ESKFState(state.p, state.v, state.q, state.ab, state.wb, P_new)
    
    # Inject error state and return
    return inject(state_updated, dx)
end

"""
    update_pos(state::ESKFState, p_meas, R_pos)

Kalman update step for a position measurement (e.g., GPS or zero position on rail).
"""
function update_pos(state::ESKFState, p_meas::SVector{3, Float64}, R_pos::SMatrix{3, 3, Float64})
    z_hat = state.p
    y = p_meas - z_hat

    H_M = @MMatrix zeros(3, 15)
    H_M[:, 1:3] = I(3)
    H = SMatrix{3, 15, Float64}(H_M)

    S = H * state.P * transpose(H) + R_pos
    K = state.P * transpose(H) * inv(S)
    dx = K * y

    I_KH = I(15) - K * H
    P_new = I_KH * state.P * transpose(I_KH) + K * R_pos * transpose(K)

    state_updated = ESKFState(state.p, state.v, state.q, state.ab, state.wb, P_new)
    return inject(state_updated, dx)
end

"""
    update_vel(state::ESKFState, v_meas, R_vel)

Kalman update step for a velocity measurement (e.g., GPS or zero velocity on rail).
"""
function update_vel(state::ESKFState, v_meas::SVector{3, Float64}, R_vel::SMatrix{3, 3, Float64})
    z_hat = state.v
    y = v_meas - z_hat

    H_M = @MMatrix zeros(3, 15)
    H_M[:, 4:6] = I(3)
    H = SMatrix{3, 15, Float64}(H_M)

    S = H * state.P * transpose(H) + R_vel
    K = state.P * transpose(H) * inv(S)
    dx = K * y

    I_KH = I(15) - K * H
    P_new = I_KH * state.P * transpose(I_KH) + K * R_vel * transpose(K)

    state_updated = ESKFState(state.p, state.v, state.q, state.ab, state.wb, P_new)
    return inject(state_updated, dx)
end

"""
    update_rate(state::ESKFState, w_meas, w_imu, R_rate)

Kalman update step for an angular rate measurement (e.g., zero angular rate on rail).
Requires the current IMU gyroscope reading (`w_imu`) to calculate the nominal predicted rate.
"""
function update_rate(state::ESKFState, w_meas::SVector{3, Float64}, w_imu::SVector{3, Float64}, R_rate::SMatrix{3, 3, Float64})
    # Nominal rate is the IMU reading minus the estimated bias
    z_hat = w_imu - state.wb
    y = w_meas - z_hat

    # The measurement z = w_true = w_imu - w_b
    # So the Jacobian w.r.t to δw_b is -I
    H_M = @MMatrix zeros(3, 15)
    H_M[:, 13:15] = -I(3)
    H = SMatrix{3, 15, Float64}(H_M)

    S = H * state.P * transpose(H) + R_rate
    K = state.P * transpose(H) * inv(S)
    dx = K * y

    I_KH = I(15) - K * H
    P_new = I_KH * state.P * transpose(I_KH) + K * R_rate * transpose(K)

    state_updated = ESKFState(state.p, state.v, state.q, state.ab, state.wb, P_new)
    return inject(state_updated, dx)
end

end