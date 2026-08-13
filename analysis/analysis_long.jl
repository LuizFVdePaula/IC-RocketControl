using RocketControl
using CairoMakie
using ControlSystems
using ISAtmosphere
using LinearAlgebra

##

rkt = stage("projects/ic_rocket_fin/ic_rocket_fin.json");
env = environment("projects/environment/test_env.json");

## Flight conditions and reference state

Vref = 100;
href = 1000;
tref = 2;
ρref = ρ_kg_m³(p_Pa(href), T_K(href));

## Static margin

trange = 0:0.1:5
Vrange = 15:1:160
static_margins = zeros(length(trange), length(Vrange))

for (i, t) in enumerate(trange)
    for (j, V) in enumerate(Vrange)
        M = V / a_m_s(T_K(href))
        XCG = RocketControl.Dynamics.calc_xcm(rkt, t) / rkt.aed.Lref
        ΔXCG = rkt.aed.XR - XCG
        static_margins[i, j] = -RocketControl.Aerodynamics.getXCP(rkt.aed, M, 1e-5, ΔXCG)
    end
end

## Plot

fig_stab = Figure(size = (600, 400), fontsize = 16)
ax_stab = Axis(fig_stab[1, 1], xlabel = "Time [s]", ylabel = "Speed [m/s]", xticks = 0:5, yticks = 0:20:160)
heatmap!(ax_stab, trange, Vrange, static_margins, colormap = :viridis)
Colorbar(fig_stab[1, 2], colorrange = extrema(static_margins))
display(fig_stab)
#save("results/static_margin.pdf", fig_stab)

##

dm = RocketControl.DynamicModel(rkt);

q_bar = 0.5 * ρref * Vref^2
m = dm.m(tref)
Jyy = dm.Jyy(tref)
Sref = dm.Sref
Lref = dm.Lref
CNα = dm.CNα
CNδq = dm.CNδq
Cmα = dm.Cmα #- ΔXCG * CNα
Cmδq = dm.Cmδq #- ΔXCG * CNδq

Mα = q_bar * Sref * Lref * Cmα / Jyy
Nα = q_bar * Sref * CNα / m
Mq = q_bar * Sref * Lref^2 * dm.Cmq / (2 * Vref * Jyy)
M_delta = q_bar * Sref * Lref * Cmδq / Jyy
N_delta = q_bar * Sref * CNδq / m

# State space model for longitudinal dynamics: x = [q, θ, γ]
A = [
    Mq   Mα        -Mα;
    1.0  0.0        0.0;
    0.0  -Nα/Vref  Nα/Vref
]

B = [M_delta; 0.0; -N_delta/Vref]

sysc = ss(A, B, [0 0 1], 0)
sysd = c2d(sysc, 0.02)

# control law given by u = -K * (x - Nx * γref) + Nu * γref
# For discrete time, steady state requires (A - I)x + Bu = 0
Maux_d = [sysd.A - I(3) sysd.B; sysd.C 0]
NxNu_d = Maux_d \ [zeros(3, 1); 1.0]
Nx_d = NxNu_d[1:3]
Nu_d = NxNu_d[4]

## Dynamic Stability: Open-Loop Roots (Discrete)
ol_poles = pole(sysd)
println("Discrete Open Loop Poles: ", ol_poles)
println("Absolute values: ", abs.(ol_poles))

## LQR Root Locus (varying R)
Q_lqr = diagm([1.3, 30, 800.0])
R_vals = 10 .^ range(1, 4, length = 300)

poles = zeros(ComplexF64, 3, length(R_vals))

R_nom = 48.0
K_nom = lqr(sysd, Q_lqr, R_nom)
cl_poles_nom = pole(ss(sysd.A - sysd.B * K_nom, sysd.B, sysd.C, sysd.D, sysd.Ts))

for (i, R_val) in enumerate(R_vals)
    K = lqr(sysd, Q_lqr, R_val)
    cl_sys = ss(sysd.A - sysd.B * K, sysd.B, sysd.C, sysd.D, sysd.Ts)
    poles[:, i] = pole(cl_sys)
end

fig_rl = Figure(size = (500, 400), fontsize = 16)
ax_rl = Axis(fig_rl[1, 1], xlabel = "Real", ylabel = "Imaginary", xticks = 0.4:0.1:1.0, yticks = -0.2:0.05:0.2, aspect = DataAspect(), limits = ((0.6, 1.05), (-0.17, 0.17)))
θ_circle = range(0, 2π, length = 200)
lines!(ax_rl, cos.(θ_circle), sin.(θ_circle), color = :gray, linestyle = :dash)
scatter!(ax_rl, real.(ol_poles), imag.(ol_poles), marker = :x, color = :black, markersize = 15, label = "Open Loop Poles")
for p in 1:3
    lines!(ax_rl, real.(poles[p, :]), imag.(poles[p, :]), linewidth = 2)
end
scatter!(ax_rl, real.(cl_poles_nom), imag.(cl_poles_nom), marker = :circle, color = :red, markersize = 15, label = "Closed Loop Poles")
r = 0.01:0.01:1
θr = sqrt(0.5) * (-log.(r))
lines!(ax_rl, r .* cos.(θr), r .* sin.(θr), color = :green, linestyle = :dot, label = "ζ = 0.707")
axislegend(ax_rl, position = :lb)
display(fig_rl)
save("results/pitch_root_locus.pdf", fig_rl)

## FPA Closed-Loop Step Response (Discrete)
# x[k+1] = (A - B*K)*x[k] + B*(K*Nx + Nu) * γref
A_cl = sysd.A - sysd.B * K_nom
B_cl = sysd.B * (K_nom * Nx_d .+ Nu_d)
sys_cl = ss(A_cl, B_cl, sysd.C, 0.0, sysd.Ts)

t_step = 0:sysd.Ts:10.0
y_step, t_out, x_step = step(sys_cl, t_step)

fig_step = Figure(size = (600, 400), fontsize = 16)
ax_step = Axis(fig_step[1, 1], xlabel = "Time [s]", ylabel = "Angle [deg]", title = "Discrete FPA Step Response")

# step() for discrete systems returns the points at exactly the sample times
lines!(ax_step, t_out, rad2deg.(y_step[:]), linewidth = 2, label = "γ (Flight Path)")
lines!(ax_step, t_out, rad2deg.(x_step[2, :]), linewidth = 2, linestyle = :dash, label = "θ (Pitch)")
lines!(ax_step, t_out, rad2deg.(ones(length(t_out))), linewidth = 1, linestyle = :dot, color = :black, label = "Reference")

axislegend(ax_step, position = :rb)
display(fig_step)
# save("results/step_response_discrete.pdf", fig_step)
