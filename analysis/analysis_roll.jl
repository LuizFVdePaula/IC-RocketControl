using RocketControl
using CairoMakie
using ControlSystems
using LinearAlgebra
using ISAtmosphere

##

rkt = stage("projects/ic_rocket_fin/ic_rocket_fin.json");
env = environment("projects/environment/test_env.json");

## Reference variables

Vref = 100;
href = 100;
tref = 2;
ρref = ρ_kg_m³(p_Pa(href), T_K(href));

dm = RocketControl.DynamicModel(rkt);

##

q_bar = 0.5 * ρref * Vref^2
m = dm.m(tref)
Jxx = dm.Jxx(tref)
Sref = dm.Sref
Lref = dm.Lref
Clp = dm.Clp
Clδp = dm.Clδp
Lp = q_bar * Sref * Lref^2 * Clp / (2 * Vref * Jxx)
Lδp = q_bar * Sref * Lref * Clδp / Jxx

Ts = 0.02
A = [0 1; 0 Lp]
B = [0; Lδp]
C = [1 0]
sysc = ss(A, B, C, 0)
sys = c2d(sysc, Ts)

##

ω_roll_desired = 2π * 1.5
ξ_roll_desired = 0.707
poles_continuous = ω_roll_desired * [-ξ_roll_desired + im * sqrt(1 - ξ_roll_desired^2), -ξ_roll_desired - im * sqrt(1 - ξ_roll_desired^2)]
poles_discrete = exp.(poles_continuous * Ts)

#continuous: Kϕ = ω_roll_desired^2 / Lδp
#continuous: Kp = (2 * ξ_roll_desired * ω_roll_desired + Lp) / Lδp

K = place(sys, poles_discrete)
K_ϕ_nom = K[1]
K_p_nom = K[2]
K_lqr = lqr(sys, diagm([1 / 3^2, 1 / 30^2]), 1 / 0.3^2)

## Root locus plot for Successive Loop Closure

fig = Figure(size = (500, 400), fontsize = 16)
ax = Axis(fig[1, 1], xlabel = "Real", ylabel = "Imaginary", aspect = DataAspect(), yticks = -0.2:0.1:0.2, limits = ((0.45, 1.05), (-0.25, 0.25)))

# Plot Unit Circle
θ_circle = range(0, 2π, length = 200)
lines!(ax, cos.(θ_circle), sin.(θ_circle), color = :gray, linestyle = :dash)

# 1. Inner Loop: Roll Rate (p) Feedback
# Plant: u -> p
sys_p = ss(sys.A, sys.B, [0 1], 0, Ts)
Kp_range = range(0, 2 * K_p_nom, length = 200)
roots_p, _, _ = rlocus(sys_p, Kp_range)

lines!(ax, real(roots_p[:, 1]), imag(roots_p[:, 1]), color = :blue, label = "Inner Loop (p)")
lines!(ax, real(roots_p[:, 2]), imag(roots_p[:, 2]), color = :blue)

# Plot the pole locations when the inner loop is closed with K_p_nom
A_inner = sys.A - sys.B * [0 K_p_nom]
poles_inner = pole(ss(A_inner, sys.B, [0 1], 0, Ts))
scatter!(ax, real(poles_inner), imag(poles_inner), marker = :diamond, color = :blue, markersize = 15, label = "Inner Closed (Kp)")

# 2. Outer Loop: Roll Angle (ϕ) Feedback
# Plant: ϕ_cmd -> ϕ, with the inner loop closed!
sys_ϕ = ss(A_inner, sys.B, [1 0], 0, Ts)
Kϕ_range = range(0, 2 * K_ϕ_nom, length=200)
roots_ϕ, _, _ = rlocus(sys_ϕ, Kϕ_range)

lines!(ax, real(roots_ϕ[:, 1]), imag(roots_ϕ[:, 1]), color = :red, label = "Outer Loop (ϕ)")
lines!(ax, real(roots_ϕ[:, 2]), imag(roots_ϕ[:, 2]), color = :red)

# Plot the final desired poles
scatter!(ax, real(poles_discrete), imag(poles_discrete), marker = :circle, color = :black, markersize = 15, label = "Final Poles (Kϕ)")

# Plot curve of ξ = 0.707
r_ξ = 0.01:0.01:1
θ_ξ = -sqrt(1 - ξ_roll_desired^2) / ξ_roll_desired * log.(r_ξ)
x_ξ = r_ξ .* cos.(θ_ξ)
y_ξ = r_ξ .* sin.(θ_ξ)
lines!(ax, x_ξ, y_ξ, label = "ξ = 0.707", linestyle = :dash, color = :green)

axislegend(ax, position = :lb)
save("results/roll_root_locus.pdf", fig)
display(fig)

