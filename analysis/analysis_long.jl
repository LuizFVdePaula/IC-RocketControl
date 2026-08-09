using RocketControl
using ISAtmosphere
using ControlSystems

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

using CairoMakie

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
sysd = c2d(sysc, 0.002)

# control law given by u = -K * (x - [0; 1; 1] * γref)

Maux = [sysc.A sysc.B; sysc.C 0]
NxNu = Maux \ [zeros(3, 1); 1.0]
Nx = NxNu[1:3]
Nu = NxNu[4]

