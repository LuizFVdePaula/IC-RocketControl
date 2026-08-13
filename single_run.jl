using RocketControl
using LinearAlgebra
using StaticArrays
using ControlSystems

## 1. Load the stage and environment
rkt = stage("projects/ic_rocket_fin/ic_rocket_fin.json");
env = environment("projects/environment/env_model.json");

# 2. Setup the ESKF noise matrices
Ts_imu = 0.002     # 500 Hz
Ts_baro = 0.010    # 100 Hz
Ts_control = 0.020 # 50 Hz control loop

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

ωn_roll = 7.0 # rad/s
ζ_roll = 0.9

target_apogee = 1000.0 # meters

# 4. Instantiate ControlParameters
params = ControlParameters(
    Ts_imu,
    Ts_baro,
    Ts_control,
    Q_eskf,
    R_baro,
    ωn_pitch,
    ζ_pitch,
    ωn_roll,
    ζ_roll,
    target_apogee
)

## Simulation

ts = range(-10.0, 20.0, step = Ts_imu) # Step matches Ts_imu exactly
ts_flight = ts[5001:end]

sim, control_hist, est, meas = simulate(rkt, env, params)
res = postprocess(rkt, env, sim, control_hist, est, meas, ts)
res_flight = res[5001:end, :]

sim_nc, control_hist_nc, est_nc, meas_nc = simulate(rkt, env, params, ctrl_func = nothing)
res_nc = postprocess(rkt, env, sim_nc, control_hist_nc, est_nc, meas_nc, ts)
res_flight_nc = res_nc[5001:end, :]

##

include("plotsim.jl")

plot_performance(res_flight, ts_flight)
plot_observer(res, ts)
plot_observer_bias(res, rkt, ts)
plot_control(res_flight, ts_flight)

##

γobs = @. atan(-res.vDobs, sqrt(res.vNobs^2 + res.vEobs^2))

fig_gamma = Figure(size = (1200, 800))

ax_gamma = Axis(fig_gamma[1, 1], xlabel = "Time [s]", ylabel = "Flight Path Angle γ [deg]")
lines!(ax_gamma, ts, rad2deg.(res.gamma), label = "True γ", linewidth = 2)
lines!(ax_gamma, ts, rad2deg.(γobs), label = "Estimated γ", linewidth = 2)
display(fig_gamma)

## Final plots

using CairoMakie

fig_final = Figure(size = (600, 900))

ax_fpa = Axis(fig_final[1, 1], ylabel = "FPA γ [deg]", limits = (nothing, (0, 90)))
lines!(ax_fpa, ts_flight, rad2deg.(res_flight.gamma), label = "Controlled")
lines!(ax_fpa, ts_flight, rad2deg.(res_flight_nc.gamma), label = "NC", linestyle = :dash)
hlines!(ax_fpa, 80.0, linestyle = :dot, color = :black, label = "Reference")
#axislegend(ax_fpa)

ax_roll = Axis(fig_final[2, 1], ylabel = "Roll Angle ϕ [deg]")
lines!(ax_roll, ts_flight, rad2deg.(res_flight.ϕ), label = "Controlled")
lines!(ax_roll, ts_flight, rad2deg.(res_flight_nc.ϕ), label = "NC", linestyle = :dash)
#axislegend(ax_roll)

ax_uq = Axis(fig_final[3, 1], ylabel = "Pitch Deflection [deg]")
lines!(ax_uq, ts_flight, rad2deg.(res_flight.δq), label = "Controlled")
lines!(ax_uq, ts_flight, rad2deg.(res_flight_nc.δq), label = "NC", linestyle = :dash)
#axislegend(ax_uq)

ax_up = Axis(fig_final[4, 1], ylabel = "Roll Deflection [deg]")
lines!(ax_up, ts_flight, rad2deg.(res_flight.δp), label = "Controlled")
lines!(ax_up, ts_flight, rad2deg.(res_flight_nc.δp), label = "NC", linestyle = :dash)
#axislegend(ax_up)

ax_V = Axis(fig_final[5, 1], ylabel = "Total Speed [m/s]")
lines!(ax_V, ts_flight, (@. sqrt(res_flight.u^2 + res_flight.v^2 + res_flight.w^2)), label = "Controlled")
lines!(ax_V, ts_flight, (@. sqrt(res_flight_nc.u^2 + res_flight_nc.v^2 + res_flight_nc.w^2)), label = "NC", linestyle = :dash)
#axislegend(ax_V)

display(fig_final)