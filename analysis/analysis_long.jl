using RocketControl
using ISAtmosphere

##

rkt = stage("projects/ic_rocket_v2/ic_rocket_v2.json");
env = environment("projects/environment/test_env.json");

## Flight conditions and reference state

Vref = 100;
href = 500;
tref = 2;
ρref = 1.1;

## Static margin

Vrange = 15:1:160
trange = 0:0.1:5
static_margins = zeros(length(Vrange), length(trange))

for (i, V) in enumerate(Vrange)
    for (j, t) in enumerate(trange)
        M = V / a_m_s(T_K(href))
        XCG = RocketControl.Dynamics.calc_xcm(rkt, t) / rkt.aed.Lref
        ΔXCG = rkt.aed.XR - XCG
        static_margins[i, j] = RocketControl.Aerodynamics.getXCP(rkt.aed, M, 1e-5, ΔXCG)
    end
end

## Plot

using CairoMakie

fig_stab = Figure(size = (600, 400), fontsize = 16)
ax_stab = Axis(fig_stab[1, 1], xlabel = "Time [s]", ylabel = "Speed [m/s]", xticks = 0:5, yticks = 0:20:160)
heatmap!(ax_stab, trange, Vrange, static_margins', colormap = :viridis)
Colorbar(fig_stab[1, 2], colorrange = extrema(static_margins))
display(fig_stab)
save("results/static_margin.pdf", fig_stab)

## Longitudinal natural frequency

