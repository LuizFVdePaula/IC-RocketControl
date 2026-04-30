using RocketControl
using ControlSystems
using LinearAlgebra

##

rkt = stage("projects/ic_rocket_v2/ic_rocket_v2.json");
env = environment("projects/environment/test_env.json");

## Flight conditions and reference state

V_cond = [90, 120];
h_cond = [100, 500];
t_cond = [1.5, 3.0]; # t₁ ∈ [0, 2], t₂ ∈ [2, 5]
ρ_cond = [1.213, 1.167]
ϕeq = deg2rad(0);
θeq = deg2rad(80);
ψeq = 0;
q0eq = RocketControl.BaseDefs.calc_q0(ϕeq, θeq, ψeq);
q1eq = RocketControl.BaseDefs.calc_q1(ϕeq, θeq, ψeq);
q2eq = RocketControl.BaseDefs.calc_q2(ϕeq, θeq, ψeq);
q3eq = RocketControl.BaseDefs.calc_q3(ϕeq, θeq, ψeq);

cond = 1
Vref = V_cond[cond]
href = h_cond[cond]
tref = t_cond[cond]
ρref = ρ_cond[cond]
xeq = [0, 0, -href, q0eq, q1eq, q2eq, q3eq, Vref, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

## Control and state estimator gain
# x = [p, q, r, α, β, δp, δq, δr]
# Linearization

#Atot = ∇(x -> dynamics(x, [0.0, 0.0, 0.0], rkt, env, tref), xeq)
#Btot = ∇(u -> dynamics(xeq, u, rkt, env, tref), [0.0, 0.0, 0.0])

model = RocketControl.GNC.DynamicModel(rkt);
sysc = RocketControl.GNC.continuousmodel(model, xeq, Vref, href, tref);
Ts = 0.01
sysd = c2d(sysc, Ts)
A = sysd.A
B = sysd.B

## LQR

Q = diagm([1 / deg2rad(1)^2, 1 / deg2rad(1)^2, 1 / deg2rad(1)^2, 1 / deg2rad(1)^2, 1 / deg2rad(1)^2]) 
Qa = diagm([1 / deg2rad(15)^2, 1 / deg2rad(15)^2, 1 / deg2rad(15)^2])
Qaug = [Q zeros(5, 3); zeros(3, 5) Qa]

R = diagm([1 / deg2rad(10)^2, 1 / deg2rad(10)^2, 1 / deg2rad(10)^2])
L = lqr(sysd, Qaug, R)

## Observer

#H = sysd.C
#M = place(A, H, [1e-1, 2e-1, 3e-2, 4e-2, 5e-2, -1e-1, -2e-1, -3e-1], :o) # Luenberger

σw = deg2rad(300) / 3
σe = deg2rad(15)
σ = [(σw)^2, σw^2, σw^2, σe^2, σe^2, σe^2] |> diagm
f = 0.5 * ρref * Vref^2 * rkt.aed.Sref
Mα = f * rkt.aed.Lref / model.Jyy(tref) * model.Cmα
Nα = f / model.m(tref) * model.CNα
Mβ = f * rkt.aed.Lref / model.Jyy(tref) * (-model.Cmα)
Nβ = f / model.m(tref) * model.CNα
G = [
    1 0         0         0 0 0
    0 Mα        0         0 0 0
    0 0         Mβ        0 0 0
    0 Nα / Vref 0         0 0 0
    0 0         Nβ / Vref 0 0 0
    0 0         0         1 0 0
    0 0         0         0 1 0
    0 0         0         0 0 1
]

#R1 = G * σ * G' * Ts
R1 = c2d(sysc, G * σ * G', Ts, opt = :o)
R2 = diagm([deg2rad(1)^2, deg2rad(1)^2, deg2rad(1)^2, 0.15^2, 0.15^2, deg2rad(1)^2, deg2rad(1)^2, deg2rad(1)^2])
K = kalman(sysd, R1, R2)

##

using CSV, DataFrames
df = CSV.read("projects/ic_rocket_v2/coefs_ic_v2_deg.csv", DataFrame)
df_rad = deepcopy(df)
df_rad.ClDELTAP = rad2deg.(df.ClDELTAP)
df_rad.CmDELTAQ = rad2deg.(df.CmDELTAQ)
df_rad.CnDELTAR = rad2deg.(df.CnDELTAR)
df_rad.CADELTAEFF2 = rad2deg.(rad2deg.(df.CADELTAEFF2))
df_rad.CYDELTAR = rad2deg.(df.CYDELTAR)
df_rad.CNDELTAQ = rad2deg.(df.CNDELTAQ)
CSV.write("projects/ic_rocket_v2/coefs_ic_v2_rad.csv", df_rad)

