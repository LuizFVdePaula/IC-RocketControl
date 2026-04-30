using RocketControl
using CairoMakie
using ControlSystems
using LinearAlgebra

##

rkt = stage("projects/ic_rocket_v2/ic_rocket_v2.json");
env = environment("projects/environment/test_env.json");

## Flight conditions and reference variables

V_cond = [90, 120];
h_cond = [100, 500];
t_cond = [1.5, 3.0];
ρ_cond = [1.213, 1.167];

cond = 1
Vref = V_cond[cond]
href = h_cond[cond]
tref = t_cond[cond]
ρref = ρ_cond[cond]

## System build
model = RocketControl.GNC.DynamicModel(rkt);
Ts = 0.01
f = 0.5 * ρref * Vref^2 * rkt.aed.Sref
Mp = f * model.Lref^2 / (2 * model.Jxx(tref) * Vref) * model.Clp
Mδp = f * model.Lref / model.Jxx(tref) * model.Clδp
τ = model.τ

A = [Mp Mδp; 0 -1 / τ]
B = [0; 1 / τ]
C = I(2)
sysc = ss(A, B, C, 0)
sys = c2d(sysc, Ts)

## LQR design

Q = diagm([1 / deg2rad(10)^2, 1 / deg2rad(1)^2])
R = 1 / deg2rad(10)^2
L = lqr(sys, Q, R)

## Kalman filter design

σ = diagm([(1.4)^2, deg2rad(10)^2])
G = I(2)
R1 = c2d(sys, G * σ * G'; opt = :o)
R2 = diagm([deg2rad(1)^2, deg2rad(1)^2])
K = kalman(sys, R1, R2; direct = true)

## Linear simulation

cont = observer_controller(sys, L, K; direct = true)
syscl = feedback(sys, cont)
resx = lsim(syscl, (x, t) -> 0, 0:0.01:0.4; x0 = [deg2rad(15), 0, 0, 0])

## Simulation result

x̂ = resx.x[3:4, :] + K * resx.x[1:2, :]

fig = Figure(size = (900, 600))

ax1 = Axis(fig[1, 1])
lines!(ax1, resx.t, resx.x[1, :], label = "p")
lines!(ax1, resx.t, x̂[1, :], label = "p̂")
axislegend(ax1)

ax2 = Axis(fig[1, 2])
lines!(ax2, resx.t, resx.x[2, :], label = "δp")
lines!(ax2, resx.t, x̂[2, :], label = "δp̂")
axislegend(ax2)

fig

## Root locus plot

#ρs = (10).^range(-4, 2, length = 100)
ρs = range(1 / deg2rad(15)^2, 1 / deg2rad(1)^2, 100)
root_locus = zeros(ComplexF64, (length(ρs), 2))

for (i, ρ) in enumerate(ρs)
    L = lqr(sys, diagm([1 / deg2rad(10)^2, ρ]), R)
    λs = eigvals(sys.A - sys.B * L)
    root_locus[i, :] = λs
end

fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1], title = "Discrete LQR Root Locus", xlabel = "Real", ylabel = "Imag", aspect = DataAspect())

θun = range(0, 2π, length = 100)
lines!(ax, cos.(θun), sin.(θun), color = :black, linestyle = :dash)

scatterlines!(ax, real(root_locus[:, 1]), imag(root_locus[:, 1]), markersize = 5, label = "λ1")
scatterlines!(ax, real(root_locus[:, 2]), imag(root_locus[:, 2]), markersize = 5, label = "λ2")

λ0 = eigvals(sys.A)
scatter!(ax, real(λ0[1]), imag(λ0[1]), marker = :circle, markersize = 10, label = "Open Loop λ₁")
scatter!(ax, real(λ0[2]), imag(λ0[2]), marker = :circle, markersize = 10, label = "Open Loop λ₂")

axislegend(ax)
fig