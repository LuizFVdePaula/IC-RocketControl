using RocketControl
using LinearAlgebra
using StaticArrays
using ControlSystems

## 1. Load the stage and environment
rkt = stage("projects/ic_rocket_fin/ic_rocket_fin.json");
env = environment("projects/environment/env_model.json");

# 2. Setup the ESKF noise matrices
Ts_imu = 0.002
Ts_baro = 0.010

# ESKF Process Noise (Q_eskf) for 15 error states:
# [δp (3), δv (3), δθ (3), δab (3), δwb (3)]
σ_p = 1e-4
σ_v = 250e-6 * 9.80665
σ_th = deg2rad(0.015)
σ_ab = 1e-5
σ_wb = 1e-5

Q_vec = [
    fill(σ_p^2, 3);
    fill(σ_v^2, 3);
    fill(σ_th^2, 3);
    fill(σ_ab^2, 3);
    fill(σ_wb^2, 3)
]
Q_eskf = SMatrix{15, 15, Float64}(diagm(Q_vec))

# ESKF Measurement Noise
R_baro = 2.0^2 # 2 m standard deviation on barometer

# 3. Setup the Attitude Autopilot Parameters
# We target a highly damped, fast response
ωn_pitch = 5.0 # rad/s
ζ_pitch = 0.7

ωn_roll = 8.0 # rad/s
ζ_roll = 0.7

target_apogee = 1000.0 # meters

# 4. Instantiate ControlParameters
params = ControlParameters(
    Ts_imu,
    Ts_baro,
    Q_eskf,
    R_baro,
    ωn_pitch,
    ζ_pitch,
    ωn_roll,
    ζ_roll,
    target_apogee
)

##

sim, control_hist, est, meas = simulate(rkt, env, params)

ts = range(-10.0, 20.0, step = Ts_imu) # Step matches Ts_imu exactly
res = postprocess(rkt, env, sim, control_hist, est, meas, ts)

# Display some final statistics
println("\n=== Simulation Results ===")
println("Max Altitude: ", maximum(res.h), " m")
println("Max Speed: ", maximum(sqrt.(res.u.^2 .+ res.v.^2 .+ res.w.^2)), " m/s")
println("Final Attitude (θ): ", rad2deg(res.θ[end]), " deg")
println("Final Estimated Attitude (θ̂): ", rad2deg(res.θobs[end]), " deg")

##

include("plotsim.jl")

plot_performance(res, ts)
plot_observer(res, ts)
plot_observer_bias(res, rkt, ts)
plot_control(res, ts)

##

γobs = @. atan(-res.vDobs, sqrt(res.vNobs^2 + res.vEobs^2))

fig_gamma = Figure(size = (1200, 800))

ax_gamma = Axis(fig_gamma[1, 1], xlabel = "Time [s]", ylabel = "Flight Path Angle γ [deg]")
lines!(ax_gamma, ts, rad2deg.(res.gamma), label = "True γ", linewidth = 2)
lines!(ax_gamma, ts, rad2deg.(γobs), label = "Estimated γ", linewidth = 2)
display(fig_gamma)