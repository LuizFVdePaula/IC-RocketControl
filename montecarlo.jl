using RocketControl
using LinearAlgebra
using StaticArrays
using ControlSystems
using CairoMakie

## 1. Load base configuration
rkt_base = stage("projects/ic_rocket_fin/ic_rocket_fin.json");
env_json = "projects/environment/env_model.json";

## 2. Setup ControlParameters
Ts_imu = 0.002
Ts_baro = 0.010
Ts_control = 0.020 # 50 Hz

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
R_baro = 2.0^2

ωn_pitch = 5.0
ζ_pitch = 0.7
ωn_roll = 5.0
ζ_roll = 0.9
target_apogee = 1000.0

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

N_runs = 100

##

println("Running $N_runs simulations WITHOUT control...")
results_uncontrolled = run_monte_carlo(N_runs, rkt_base, env_json, params; ctrl_func=nothing)

println("Running $N_runs simulations WITH FPA control...")
results_controlled = run_monte_carlo(N_runs, rkt_base, env_json, params; ctrl_func=control_fpa)

## 3. Post-process to extract key metrics
# Final East (y), North (x) positions for dispersion
pos_unc = [(df.y[end], df.x[end]) for df in results_uncontrolled]
pos_ctrl = [(df.y[end], df.x[end]) for df in results_controlled]

apogee_unc = [maximum(df.h) for df in results_uncontrolled]
apogee_ctrl = [maximum(df.h) for df in results_controlled]

# 4. Plotting
fig = Figure(size = (1200, 600), fontsize = 18)

# Plot 1: Apogee Distribution
ax_apogee = Axis(fig[1, 1], xlabel = "Apogee [m]", ylabel = "Count", title = "Apogee Dispersion")
hist!(ax_apogee, apogee_unc, bins = 15, color = (:red, 0.5), label = "Uncontrolled")
hist!(ax_apogee, apogee_ctrl, bins = 15, color = (:blue, 0.5), label = "Controlled (FPA)")
axislegend(ax_apogee)

# Plot 2: Landing Dispersion (East vs North)
ax_disp = Axis(fig[1, 2], xlabel = "East [m]", ylabel = "North [m]", title = "Landing Dispersion", aspect = DataAspect())
scatter!(ax_disp, [p[1] for p in pos_unc], [p[2] for p in pos_unc], color = (:red, 0.7), label = "Uncontrolled", markersize = 8)
scatter!(ax_disp, [p[1] for p in pos_ctrl], [p[2] for p in pos_ctrl], color = (:blue, 0.7), label = "Controlled (FPA)", markersize = 8)
axislegend(ax_disp, position = :rt)

display(fig)

# Save the plot
save("results/montecarlo_dispersion.pdf", fig)
println("Monte Carlo analysis complete! Plots saved to results/montecarlo_dispersion.pdf")

##

fig_alts = Figure(size = (600, 400))
ax_alts = Axis(fig_alts[1, 1], xlabel = "Time [s]", ylabel = "FPA γ [deg]")
ax_rolls = Axis(fig_alts[2, 1], xlabel = "Time [s]", ylabel = "Roll Angle ϕ [deg]")
for result ∈ results_uncontrolled
    idx = argmax(result.h)
    lines!(ax_alts, ts[5001:idx], rad2deg.(result.gamma[5001:idx]), color = :red, alpha = 0.05)
    lines!(ax_rolls, ts[5001:idx], rad2deg.(result.ϕ[5001:idx]), color = :red, alpha = 0.05)
end
for result ∈ results_controlled
    idx = argmax(result.h)
    lines!(ax_alts, ts[5001:idx], rad2deg.(result.gamma[5001:idx]), color = :blue, alpha = 0.05)
    lines!(ax_rolls, ts[5001:idx], rad2deg.(result.ϕ[5001:idx]), color = :blue, alpha = 0.05)
end
display(fig_alts)