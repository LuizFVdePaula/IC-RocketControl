module Propulsion

export PropulsionSubsystem, SolidPropellant, SolidPropulsion,
thrust, propellant_mass, propellant_mass_derivative, propellant_center_of_mass,
propellant_inertia_tensor, propellant_inertia_tensor_derivative

using CSV
using DataFrames
using Interpolations
using LinearAlgebra
using StaticArrays

abstract type PropulsionSubsystem end

"""
    function thrust(prp::PropulsionSubsystem, t)

Obtain current thrust at time `t`.
"""
function thrust end

function propellant_mass end

function propellant_mass_derivative end

function propellant_center_of_mass end

function propellant_inertia_tensor end

function propellant_inertia_tensor_derivative end

struct SolidPropellant{T}
    m0::Float64
    ρ::Float64
    L::Float64
    R::Float64
    m::T
    mdot::T
    Ixx::T
    Iyy::T
    Ixxdot::T
    Iyydot::T
end

function SolidPropellant(T::AbstractVector, t::AbstractVector, It, m0, ρ, L, R)
    mdot = -T * m0 / It
    m = m0 .+ [0; cumsum((mdot[begin:end-1] + mdot[begin+1:end]) .* diff(t))] / 2
    Ixx = m * R^2 - m.^2 / (2π * ρ * L)
    Iyy = m * (R^2 / 2 + L^2 / 12) - m.^2 / (4π * ρ * L)
    Ixxdot = mdot * R^2 - m .* mdot / (π * ρ * L)
    Iyydot = mdot * (R^2 / 2 + L^2 / 12) - m .* mdot / (2π * ρ * L)
    return SolidPropellant(
        m0, ρ, L, R,
        extrapolate(interpolate((t,), m, Gridded(Linear())), Flat()),
        extrapolate(interpolate((t,), mdot, Gridded(Linear())), Flat()),
        extrapolate(interpolate((t,), Ixx, Gridded(Linear())), Flat()),
        extrapolate(interpolate((t,), Iyy, Gridded(Linear())), Flat()),
        extrapolate(interpolate((t,), Ixxdot, Gridded(Linear())), Flat()),
        extrapolate(interpolate((t,), Iyydot, Gridded(Linear())), Flat())
    )
end

struct SolidPropulsion{P, S} <: PropulsionSubsystem
    propellant::P
    T::S
    Pe::Float64
    It::Float64
    tb::Float64
    Ae::Float64
    xcm::Float64
    re::SVector{3, Float64}
    ξ::Float64
    η::Float64
end

function SolidPropulsion(T::AbstractVector, t::AbstractVector, m0, ρ, L, R, Pe, Ae, xcm, re, ξ, η)
    It = sum((T[begin:end-1] + T[begin+1:end]) .* diff(t)) / 2
    propellant = SolidPropellant(T, t, It, m0, ρ, L, R)
    return SolidPropulsion(
        propellant,
        extrapolate(interpolate((t,), T, Gridded(Linear())), 0.0),
        Pe, It, last(t), Ae, xcm, re, ξ, η
    )
end

function thrust(prp::SolidPropulsion, t)
    if !(0.0 <= t <= prp.tb)
        return 0.0
    end
    return prp.T(t)
end

function propellant_mass(prp::SolidPropulsion, t)
    if t < 0.0
        return prp.propellant.m0
    end
    prp.propellant.m(t)
end

function propellant_mass_derivative(prp::SolidPropulsion, t)
    if !(0.0 <= t <= prp.tb)
        return 0.0
    end
    return prp.propellant.mdot(t)
end

propellant_center_of_mass(prp::SolidPropulsion, t) = prp.xcm

function propellant_inertia_tensor(prp::SolidPropulsion, t)
    Ixx = prp.propellant.Ixx(t)
    Iyy = prp.propellant.Iyy(t)
    J₀ = SMatrix{3, 3}(diagm([Ixx, Iyy, Iyy]))
    return J₀
end

function propellant_inertia_tensor_derivative(prp::SolidPropulsion, t)
    if !(0.0 <= t <= prp.tb)
        return zeros(SMatrix{3, 3})
    end
    İxx = prp.propellant.Ixxdot(t)
    İyy = prp.propellant.Iyydot(t)
    J̇₀ = SMatrix{3, 3}(diagm([İxx, İyy, İyy]))
    return J̇₀
end

function from_dict(dict::AbstractDict)
    filepath = dict["thrust_curve"]
    m0 = dict["propellant_mass"]
    Pe = dict["exit_pressure"]
    Ae = dict["exit_area"]
    xcm = dict["center_of_mass"]
    xe = dict["exit_position_x"]
    ye = dict["exit_position_y"]
    ze = dict["exit_position_z"]
    ρ = dict["propellant_density"]
    L = dict["propellant_length"]
    R = dict["propellant_diameter"] / 2
    ξ = deg2rad(dict["deflection_angle_z"])
    η = deg2rad(dict["deflection_angle_y"])
    re = SVector(xe, ye, ze)
    curve = CSV.read(filepath, DataFrame)
    return SolidPropulsion(curve.thrust, curve.time, m0, ρ, L, R, Pe, Ae, xcm, re, ξ, η)
end

function from_montecarlo(prp::SolidPropulsion, σT, σreyz, σξη)
    fT = 1 + randn() * σT
    T_mc = prp.T.itp.coefs .* fT
    reyz = randn() * σreyz
    θreyz = rand() * 2π
    ξη = randn() * σξη
    θξη = rand() * 2π
    rey, rez = reyz .* sincos(θreyz)
    re = SVector(prp.re[1], rey, rez)
    ξ, η = ξη .* sincos(θξη)
    return SolidPropulsion(
        T_mc,
        prp.T.itp.knots[1],
        prp.propellant.m0,
        prp.propellant.ρ,
        prp.propellant.L,
        prp.propellant.R,
        prp.Pe,
        prp.Ae,
        prp.xcm,
        re,
        ξ,
        η
    )
end

end