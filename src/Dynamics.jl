module Dynamics

export aeroloads, proploads, loads, dynamics
export calc_mdot, calc_xcm, calc_J, calc_Jdot, calc_ρ, calc_mach, calc_αT, calc_ϕA

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
    J_str = stg.str.J + stg.str.m * (stg.str.xcm - xcm)^2 * Myz # review
    J_prp = J₀_prp + m_prp * (xcm_prp - xcm)^2 * Myz
    J = J_str + J_prp
    return J
end

function calc_Jdot(stg::Stage, t_stg, xcm = calc_xcm(stg, t_stg))
    Myz = SMatrix{3, 3}([0 0 0; 0 1 0; 0 0 1])
    J̇₀_prp = propellant_inertia_tensor_derivative(stg.prp, t_stg)
    ṁ_prp = propellant_mass_derivative(stg.prp, t_stg)
    xcm_prp = propellant_center_of_mass(stg.prp, t_stg)
    J̇ = J̇₀_prp + ṁ_prp * (xcm - xcm_prp)^2 * Myz # review
    return J̇
end

calc_ρ(h) = ρ_kg_m³(p_Pa(h), T_K(h))

calc_mach(vBA, h) = abs(vBA) / a_m_s(T_K(h))

calc_αT(vBA) = atan(norm(vBA[2:3]) / vBA[1])

calc_ϕA(vBA) = mod2pi(atan(vBA[2], vBA[3]))

function aeroloads(h, vBA, ωBA, u, stg, t_stg)
    vBAnorm = norm(vBA)
    q̄ = 0.5 * calc_ρ(h) * vBAnorm^2
    M = calc_mach(vBAnorm, h)
    αT = calc_αT(vBA)
    ϕA = calc_ϕA(vBA)
    TBR = rotX(ϕA)
    δ = transpose(TBR) * u

    S = stg.aed.Sref
    CA = t_stg < stg.prp.tb ? getCAon(stg.aed, M, αT, ϕA, δ[2], δ[3]) : getCAoff(stg.aed, M, αT, ϕA, δ[2], δ[3])
    CY = getCY(stg.aed, M, αT, ϕA, δ[3])
    CN = getCN(stg.aed, M, αT, ϕA, δ[2])
    Faero = q̄ * S * TBR * SVector(CA, CY, CN)

    L = stg.aed.Lref
    XCG = calc_xcm(stg, t_stg) / L
    ΔXCG = stg.aed.XR - XCG
    Ωnorm = transpose(TBR) * ωBA * L / (2 * vBAnorm)
    Cl = getCl(stg.aed, M, αT, ϕA, Ωnorm[1], δ[1])
    Cm = getCm(stg.aed, M, αT, ϕA, ΔXCG, Ωnorm[2], δ[2])
    Cn = getCn(stg.aed, M, αT, ϕA, ΔXCG, Ωnorm[3], δ[3])
    Maero = q̄ * S * L * TBR * SVector(Cl, Cm, Cn)

    return (Faero, Maero)
end

function proploads(h, stg, t_stg, xcm = calc_xcm(stg, t_stg))
    T = thrust(stg, t_stg)
    ΔP = stg.Pe - p_Pa(h)
    Ae = stg.Ae
    ξ = stg.ξ
    η = stg.η
    Fprop = (T + ΔP * Ae) * SVector(cos(η) * cos(ξ), cos(η) * sin(ξ), -sin(η))
    re = stg.re - SVector(xcm, 0, 0)
    Mprop = re × Fprop
    return (Fprop, Mprop)
end

function loads(sv, u, stg::Stage, t_stg)
    h = -sv[3]
    vBG = SVector(sv[8:10])
    ωBG = SVector(sv[11:13])
    vBA = vBG #- TBG * windspeed() # wrap wind speed in interface
    ωBA = ωBG # wind does not rotate... or does it?
    (Faero, Maero) = aeroloads(h, vBA, ωBA, u, stg, t_stg)
    (Fprop, Mprop) = proploads(h, stg, t_stg)
    F = Faero + Fprop
    M = Maero + Mprop
    return (F, M)
end

"""
    function dynamics
"""
function dynamics(sv, u, stg, env, t)
    #sv: [x, y, z, q0, q1, q2, q3, u, v, w, p, q, r]
    #u: [δp, δq, δr]

    TBG = rotXYZ(sv[4], sv[5], sv[6], sv[7])
    vBG = SVector(sv[8:10])
    ωBG = SVector(sv[11:13])

    # kinematic equations
    xyzdot = transpose(TBG) * vBG
    p, q, r = ωBG
    Ωquat = SMatrix{4, 4}([0 -p -q -r; p 0 r -q; q -r 0 p; r q -p 0])
    quat = SVector(sv[4:7])
    quatdot = 0.5 * Ωquat * quat - 0.5 * quat * (1 - 1 / (transpose(quat) * quat))

    # dynamic equations
    (F, M) = loads(sv, u, stg, t)
    g = TBG * SVector(0, 0, env.g)
    m = stage_mass(stg, t)
    ṁ = calc_mdot(stg.prp, t)
    xcm = calc_xcm(stg, t)
    re = stg.prp.re - SVector(xcm, 0, 0)
    J = calc_J(stg, t, xcm)
    J̇ = calc_Jdot(stg, t, xcm)
    uvwdot = -ωBG × vBG + F / m + g
    pqrdot = J \ (-ωBG × (J * ωBG) + M - J̇ * ωBG - ṁ * re × (ωBG × re))

    return SVector([xyzdot; quatdot; uvwdot; pqrdot])
end

end