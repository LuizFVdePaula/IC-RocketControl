module Dynamics

export aeroloads, proploads, loads, dynamics, ∇
export calc_mdot, calc_xcm, calc_J, calc_Jdot, calc_ρ, calc_mach, calc_αT, calc_ϕA, stage_mass

using ..BaseDefs, ..EnvironmentDefs, ..StageDefs, ..Aerodynamics, ..Propulsion
using ISAtmosphere
using LinearAlgebra
using StaticArrays

function stage_mass(stg::Stage, t_stg)
    return stg.str.m + propellant_mass(stg.prp, t_stg)
end

function calc_mdot(prp::SolidPropulsion, t_stg)
    return propellant_mass_derivative(prp, t_stg)
end

function calc_xcm(stg::Stage, t_stg)
    m_str = stg.str.m
    xcm_str = stg.str.xcm
    m_prp = propellant_mass(stg.prp, t_stg)
    xcm_prp = propellant_center_of_mass(stg.prp, t_stg)
    return (m_str * xcm_str + m_prp * xcm_prp) / (m_str + m_prp)
end

function calc_J(stg::Stage, t_stg, xcm = calc_xcm(stg, t_stg))
    Myz = SMatrix{3, 3}([0 0 0; 0 1 0; 0 0 1])
    J₀_prp = propellant_inertia_tensor(stg.prp, t_stg)
    m_prp = propellant_mass(stg.prp, t_stg)
    xcm_prp = propellant_center_of_mass(stg.prp, t_stg)
    J_str = stg.str.J + stg.str.m * (stg.str.xcm - xcm)^2 * Myz
    J_prp = J₀_prp + m_prp * (xcm_prp - xcm)^2 * Myz
    J = J_str + J_prp
    return J
end

function calc_Jdot(stg::Stage, t_stg, xcm = calc_xcm(stg, t_stg))
    Myz = SMatrix{3, 3}([0 0 0; 0 1 0; 0 0 1])
    J̇₀_prp = propellant_inertia_tensor_derivative(stg.prp, t_stg)
    ṁ_prp = propellant_mass_derivative(stg.prp, t_stg)
    xcm_prp = propellant_center_of_mass(stg.prp, t_stg)
    J̇ = J̇₀_prp + ṁ_prp * (xcm - xcm_prp)^2 * Myz
    return J̇
end

calc_ρ(h) = ρ_kg_m³(p_Pa(h), T_K(h))

calc_mach(vBA, h) = abs(vBA) / a_m_s(T_K(h))

calc_αT(vBA) = atan(norm(vBA[2:3]) / vBA[1])

calc_ϕA(vBA) = mod2pi(atan(vBA[2], vBA[3]))

function aeroloads(h, vBA, ωBA, δ, stg, t_stg)
    vBAnorm = norm(vBA)
    q̄ = 0.5 * calc_ρ(h) * vBAnorm^2
    M = calc_mach(vBAnorm, h)
    αT = calc_αT(vBA)
    ϕA = calc_ϕA(vBA)
    L = stg.aed.Lref
    S = stg.aed.Sref
    XCG = calc_xcm(stg, t_stg) / L
    ΔXCG = stg.aed.XR - XCG
    ωBAnorm = ωBA * L / (2 * vBAnorm)
    on = t_stg < stg.prp.tb

    Fcoef, Mcoef = aerodynamic_coefficients(stg.aed, M, αT, ϕA, ΔXCG, ωBAnorm, δ, on)
    Faero = q̄ * S * Fcoef
    Maero = q̄ * S * L * Mcoef

    return (Faero, Maero)
end

function proploads(h, stg, t_stg, xcm = calc_xcm(stg, t_stg))
    T = thrust(stg.prp, t_stg)
    ΔP = stg.prp.Pe - p_Pa(h)
    Ae = stg.prp.Ae
    ξ = stg.prp.ξ
    η = stg.prp.η
    Fprop = (T + ΔP * Ae) * SVector(cos(η) * cos(ξ), cos(η) * sin(ξ), -sin(η))
    re = stg.prp.re - SVector(xcm, 0, 0)
    Mprop = re × Fprop
    return (Fprop, Mprop)
end

function loads(sv, δ, stg::Stage, env, t_stg, TBG)
    h = -sv[3]
    vBG = SVector{3}(sv[8:10])
    ωBG = SVector{3}(sv[11:13])
    vBA = vBG - TBG * windspeed(env, h)
    ωBA = ωBG # wind does not rotate... or does it?
    (Faero, Maero) = aeroloads(h, vBA, ωBA, δ, stg, t_stg)
    (Fprop, Mprop) = proploads(h, stg, t_stg)
    F = Faero + Fprop
    M = Maero + Mprop
    return (F, M)
end

function ∇(f::Function, x₀::AbstractVector)
    m = length(f(x₀))
    n = length(x₀)
    F = zeros(m, n)
    ε = 1e-5
    for col in 1:n
        x₊ = x₀ |> copy |> float
        x₊[col] += ε
        x₋ = x₀ |> copy |> float
        x₋[col] -= ε
        F[:, col] = (f(x₊) - f(x₋)) / 2ε
    end
    return F
end

function ∇(f::Function, x₀::Real)
    m = length(f(x₀))
    F = zeros(m)
    ε = 1e-5
    x₊ = x₀ |> copy |> float
    x₊ += ε
    x₋ = x₀ |> copy |> float
    x₋ -= ε
    F = (f(x₊) - f(x₋)) / 2ε
    return F
end

"""
    dynamics(sv, δ, stg, env, t)

Obtain the time derivative of state vector `sv` subject to deflections `δ` at time instant `t`.
"""
function dynamics(sv, u, stg, env, t)
    #sv: [x, y, z, q0, q1, q2, q3, u, v, w, p, q, r, δp, δq, δr, δ̇p, δ̇q, δ̇r]
    #u: [up, uq, ur]

    TBG = rotXYZ(sv[4], sv[5], sv[6], sv[7])
    vBG = SVector{3}(sv[8:10])
    ωBG = SVector{3}(sv[11:13])

    # kinematic equations
    xyzdot = transpose(TBG) * vBG
    p, q, r = ωBG
    Ωquat = SMatrix{4, 4}([0 -p -q -r; p 0 r -q; q -r 0 p; r q -p 0])
    quat = SVector{4}(sv[4:7])
    quatdot = 0.5 * Ωquat * quat - 0.5 * quat * (1 - 1 / (transpose(quat) * quat))

    # dynamic equations
    δ = SVector{3}(sv[14:16])
    (F, M) = loads(sv, δ, stg, env, t, TBG)
    g = TBG * SVector(0, 0, gravity)
    m = stage_mass(stg, t)
    ṁ = calc_mdot(stg.prp, t)
    xcm = calc_xcm(stg, t)
    re = stg.prp.re - SVector(xcm, 0, 0)
    J = calc_J(stg, t, xcm)
    J̇ = calc_Jdot(stg, t, xcm)
    uvwdot = -ωBG × vBG + F / m + g
    pqrdot = J \ (-ωBG × (J * ωBG) + M - J̇ * ωBG + ṁ * re × (ωBG × re))
    ωn = 70.0 # TODO insert as system input
    ξ  = 1.0  # TODO insert as system input
    δdot = sv[17:19]
    δdotdot = -2 * ξ * ωn * δdot + ωn^2 * (SVector{3}(u) - δ)
    # rail constraints
    if -sv[3] < 5.0
        quatdot = SVector(0, 0, 0, 0)
        uvwdot = SVector(max(uvwdot[1], 0), 0, 0)
        pqrdot = SVector(0, 0, 0)
    end

    return SVector{19}([xyzdot; quatdot; uvwdot; pqrdot; δdot; δdotdot])
end

end