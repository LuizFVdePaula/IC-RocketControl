module Aerodynamics

export ActiveAerodynamics
export getCAon, getCAoff, getCY, getCN, getCl, getCm, getCn, getXCP

using CSV, DataFrames, Interpolations

"""
    abstract type BaseCoefficientsModel

A base model representing aerodynamic coefficients without deflection (δ = 0).
"""
abstract type BaseCoefficientsModel end

"""
    struct InterpolatedBaseModel

Aerodynamic base coefficients represented as interpolations of (`M`, `αT`, `ϕA`).

The coefficients are given in aeroballistic coordinate system.
"""
struct InterpolatedBaseModel{T} <: BaseCoefficientsModel
    CAon::T
    CAoff::T
    CY::T
    CN::T
    Cl::T
    Cm::T
    Cn::T
    Clp::T
    Cmq::T
    Cnr::T
end

function InterpolatedBaseModel(M, αT, ϕA, coeffs, scheme)
    CAon = interp_coeff(M, αT, ϕA, coeffs.CAon, scheme)
    CAoff = interp_coeff(M, αT, ϕA, coeffs.CAoff, scheme)
    CY = interp_coeff(M, αT, ϕA, coeffs.CY, scheme)
    CN = interp_coeff(M, αT, ϕA, coeffs.CN, scheme)
    Cl = interp_coeff(M, αT, ϕA, coeffs.Cl, scheme)
    Cm = interp_coeff(M, αT, ϕA, coeffs.Cm, scheme)
    Cn = interp_coeff(M, αT, ϕA, coeffs.Cn, scheme)
    Clp = interp_coeff(M, αT, ϕA, coeffs.ClP, scheme)
    Cmq = interp_coeff(M, αT, ϕA, coeffs.CmQ, scheme)
    Cnr = interp_coeff(M, αT, ϕA, coeffs.CnR, scheme)
    return InterpolatedBaseModel(CAon, CAoff, CY, CN, Cl, Cm, Cn, Clp, Cmq, Cnr)
end

getCAon(b::InterpolatedBaseModel, M, αT, ϕA) = b.CAon(M, αT, ϕA)

getCAoff(b::InterpolatedBaseModel, M, αT, ϕA) = b.CAoff(M, αT, ϕA)

getCY(b::InterpolatedBaseModel, M, αT, ϕA) = b.CY(M, αT, ϕA)

getCN(b::InterpolatedBaseModel, M, αT, ϕA) = b.CN(M, αT, ϕA)

getCl(b::InterpolatedBaseModel, M, αT, ϕA, p) = b.Cl(M, αT, ϕA) + b.Clp(M, αT, ϕA) * p

getCm(b::InterpolatedBaseModel, M, αT, ϕA, ΔXCG, q) = b.Cm(M, αT, ϕA) + b.Cmq(M, αT, ϕA) * q - b.CN(M, αT, ϕA) * ΔXCG

getCn(b::InterpolatedBaseModel, M, αT, ϕA, ΔXCG, r) = b.Cn(M, αT, ϕA) + b.Cnr(M, αT, ϕA) * r + b.CY(M, αT, ϕA) * ΔXCG

getXCP(b::InterpolatedBaseModel, M, αT, ϕA, ΔXCG) = -getCm(b, M, αT, ϕA, ΔXCG, 0) / getCN(b, M, αT, ϕA)

"""
    struct SymmetricBaseModel

Aerodynamic base coefficients represented as interpolations of (`M`, `αT`) and periodic functions of `ϕA`.

The coefficients are given in aeroballistic coordinate system.
"""
struct SymmetricBaseModel{T} <: BaseCoefficientsModel
    CAon::T
    CAoff::T
    CYϕ::T
    CN0::T
    CNϕ::T
    Clϕ::T
    Cm0::T
    Cmϕ::T
    Cnϕ::T
    Clp::T
    Cmq::T
    Cnr::T
end

function SymmetricBaseModel(M, αT, coeff, scheme)
    CAon = interp_coeff(M, αT, coeff.CAon, scheme)
    CAoff = interp_coeff(M, αT, coeff.CAoff, scheme)
    CYϕ = interp_coeff(M, αT, coeff.CYPHI, scheme)
    CN0 = interp_coeff(M, αT, coeff.CN0, scheme)
    CNϕ = interp_coeff(M, αT, coeff.CNPHI, scheme)
    Clϕ = interp_coeff(M, αT, coeff.ClPHI, scheme)
    Cm0 = interp_coeff(M, αT, coeff.Cm0, scheme)
    Cmϕ = interp_coeff(M, αT, coeff.CmPHI, scheme)
    Cnϕ = interp_coeff(M, αT, coeff.CnPHI, scheme)
    Clp = interp_coeff(M, αT, coeff.ClP, scheme)
    Cmq = interp_coeff(M, αT, coeff.CmQ, scheme)
    Cnr = interp_coeff(M, αT, coeff.CnR, scheme)
    return SymmetricBaseModel(CAon, CAoff, CYϕ, CN0, CNϕ, Clϕ, Cm0, Cmϕ, Cnϕ, Clp, Cmq, Cnr)
end

getCAon(b::SymmetricBaseModel, M, αT, _) = b.CAon(M, αT)

getCAoff(b::SymmetricBaseModel, M, αT, _) = b.CAoff(M, αT)

getCY(b::SymmetricBaseModel, M, αT, ϕA) = b.CYϕ(M, αT) * sin(4ϕA)

getCN(b::SymmetricBaseModel, M, αT, ϕA) = b.CN0(M, αT) + b.CNϕ(M, αT) * sin(2ϕA)^2

getCl(b::SymmetricBaseModel, M, αT, ϕA, p) = b.Clϕ(M, αT) * sin(4ϕA) + b.Clp(M, αT) * p

getCm(b::SymmetricBaseModel, M, αT, ϕA, ΔXCG, q) = b.Cm0(M, αT) + b.Cmϕ(M, αT) * sin(2ϕA)^2 + b.Cmq(M, αT) * q - getCN(b, M, αT, ϕA) * ΔXCG

getCn(b::SymmetricBaseModel, M, αT, ϕA, ΔXCG, r) = b.Cnϕ(M, αT) * sin(4ϕA) + b.Cnr(M, αT) * r + getCY(b, M, αT, ϕA) * ΔXCG

getXCP(b::SymmetricBaseModel, M, αT, ϕA, ΔXCG) = -getCm(b, M, αT, ϕA, ΔXCG, 0) / getCN(b, M, αT, ϕA)

abstract type DeflectionCoefficientsModel end

"""
    struct InterpolatedDeflectionModel

Aerodynamic control coefficients represented as interpolations of (`M`, `αT`, `ϕA`).

The coefficients are given in aeroballistic coordinate system.
"""
struct InterpolatedDeflectionModel{T} <: DeflectionCoefficientsModel
    CAδeff2::T
    CYδr::T
    CNδq::T
    Clδp::T
    Cmδq::T
    Cnδr::T
end

function InterpolatedDeflectionModel(M, αT, ϕA, coeffs, scheme)
    CAδeff2 = interp_coeff(M, αT, ϕA, coeffs.CADELTAEFF2, scheme)
    CYδr = interp_coeff(M, αT, ϕA, coeffs.CYDELTAR, scheme)
    CNδq = interp_coeff(M, αT, ϕA, coeffs.CNDELTAQ, scheme)
    Clδp = interp_coeff(M, αT, ϕA, coeffs.ClDELTAP, scheme)
    Cmδq = interp_coeff(M, αT, ϕA, coeffs.CmDELTAQ, scheme)
    Cnδr = interp_coeff(M, αT, ϕA, coeffs.CnDELTAR, scheme)
    return InterpolatedDeflectionModel(CAδeff2, CYδr, CNδq, Clδp, Cmδq, Cnδr)
end

getCAδeff2(d::InterpolatedDeflectionModel, M, αT, ϕA) = d.CAδeff2(M, αT, ϕA)

getCYδr(d::InterpolatedDeflectionModel, M, αT, ϕA) = d.CYδr(M, αT, ϕA)

getCNδq(d::InterpolatedDeflectionModel, M, αT, ϕA) = d.CNδq(M, αT, ϕA)

getClδp(d::InterpolatedDeflectionModel, M, αT, ϕA) = d.Clδp(M, αT, ϕA)

getCmδq(d::InterpolatedDeflectionModel, M, αT, ϕA) = d.Cmδq(M, αT, ϕA)

getCnδr(d::InterpolatedDeflectionModel, M, αT, ϕA) = d.Cnδr(M, αT, ϕA)

struct SymmetricDeflectionModel{T} <: DeflectionCoefficientsModel
    CAδeff2::T
    CYδr::T
    CNδq::T
    Clδp::T
    Cmδq::T
    Cnδr::T
end

function SymmetricDeflectionModel(M, αT, coeffs, scheme)
    CAδeff2 = interp_coeff(M, αT, coeffs.CADELTAEFF2, scheme)
    CYδr = interp_coeff(M, αT, coeffs.CYDELTAR, scheme)
    CNδq = interp_coeff(M, αT, coeffs.CNDELTAQ, scheme)
    Clδp = interp_coeff(M, αT, coeffs.ClDELTAP, scheme)
    Cmδq = interp_coeff(M, αT, coeffs.CmDELTAQ, scheme)
    Cnδr = interp_coeff(M, αT, coeffs.CnDELTAR, scheme)
    return SymmetricDeflectionModel(CAδeff2, CYδr, CNδq, Clδp, Cmδq, Cnδr)
end

getCAδeff2(d::SymmetricDeflectionModel, M, αT, _) = d.CAδeff2(M, αT)

getCYδr(d::SymmetricDeflectionModel, M, αT, _) = d.CYδr(M, αT)

getCNδq(d::SymmetricDeflectionModel, M, αT, _) = d.CNδq(M, αT)

getClδp(d::SymmetricDeflectionModel, M, αT, _) = d.Clδp(M, αT)

getCmδq(d::SymmetricDeflectionModel, M, αT, _) = d.Cmδq(M, αT)

getCnδr(d::SymmetricDeflectionModel, M, αT, _) = d.Cnδr(M, αT)

struct ActiveAerodynamics{B<:BaseCoefficientsModel, D<:DeflectionCoefficientsModel}
    Lref::Float64
    Sref::Float64
    XR::Float64
    base::B
    deflection::D
end

getCAon(aed::ActiveAerodynamics, M, αT, ϕA, δq, δr) = getCAon(aed.base, M, αT, ϕA) + getCAδeff2(aed.deflection, M, αT, ϕA) * (abs(δq) + abs(δr))^2 / 4

getCAoff(aed::ActiveAerodynamics, M, αT, ϕA, δq, δr) = getCAoff(aed.base, M, αT, ϕA) + getCAδeff2(aed.deflection, M, αT, ϕA) * (abs(δq) + abs(δr))^2 / 4

getCY(aed::ActiveAerodynamics, M, αT, ϕA, δr) = getCY(aed.base, M, αT, ϕA) + getCYδr(aed.deflection, M, αT, ϕA) * δr

getCN(aed::ActiveAerodynamics, M, αT, ϕA, δq) = getCN(aed.base, M, αT, ϕA) + getCNδq(aed.deflection, M, αT, ϕA) * δq

getCl(aed::ActiveAerodynamics, M, αT, ϕA, p, δp) = getCl(aed.base, M, αT, ϕA, p) + getClδp(aed.deflection, M, αT, ϕA) * δp

"""
    function getCm(aed::ActiveAerodynamics, M, αT, ϕA, ΔXCG, q, δq)
    
Obtain pitch moment coefficient at condition of `M`, `αT`, `ϕA`, `ΔXCG`, `q`, `δq`.

Description of inputs:
- `M`: mach number.
- `αT`: total angle of attack.
- `ϕA`: aerodynamic roll angle.
- `ΔXCG`: normalized position difference between reference point and center of gravity (ΔXCG = XR - XCG).
- `q`: normalized pitch rate in aeroballistic coordinate system.
- `δq`: pitch deflection in aeroballistic coordinate system.
"""
getCm(aed::ActiveAerodynamics, M, αT, ϕA, ΔXCG, q, δq) = getCm(aed.base, M, αT, ϕA, ΔXCG, q) + getCmδq(aed.deflection, M, αT, ϕA) * δq

"""
    function getCn(aed::ActiveAerodynamics, M, αT, ϕA, ΔXCG, r, δr)
    
Obtain yaw moment coefficient at condition of `M`, `αT`, `ϕA`, `ΔXCG`, `r`, `δr`.

Description of inputs:
- `M`: mach number.
- `αT`: total angle of attack.
- `ϕA`: aerodynamic roll angle.
- `ΔXCG`: normalized position difference between reference point and center of gravity (ΔXCG = XR - XCG).
- `r`: normalized yaw rate in aeroballistic coordinate system.
- `δr`: yaw deflection in aeroballistic coordinate system.
"""
getCn(aed::ActiveAerodynamics, M, αT, ϕA, ΔXCG, r, δr) = getCn(aed.base, M, αT, ϕA, ΔXCG, r) + getCnδr(aed.deflection, M, αT, ϕA) * δr

"""
    function getXCP(aed::ActiveAerodynamics, M, αT, ϕA, ΔXCG)

Obtain normalized center of pressure w.r.t. center of mass at condition of `M`, `αT`, `ϕA`, `ΔXCG`.

Description of inputs:
- `M`: mach number.
- `αT`: total angle of attack.
- `ϕA`: aerodynamic roll angle.
- `ΔXCG`: normalized position difference between reference point and center of gravity (ΔXCG = XR - XCG).

Since the value is given normalized w.r.t. center of mass, the static margin is given in calibers by SM = -XCP.
"""
getXCP(aed::ActiveAerodynamics, M, αT, ϕA, ΔXCG) = getXCP(aed.base, M, αT, ϕA, ΔXCG)

"""
    from_dict(dict::AbstractDict)

Construct `ActiveAerodynamics` from `dict` containing reference length and area, and path to aerodynamic coefficients.

Useful when reading `.json` files.
"""
function from_dict(dict::AbstractDict)
    Lref = dict["reference_length"]
    Sref = dict["reference_area"]
    XR = dict["reference_position"] / Lref
    model = dict["model"]
    extrapolation = dict["extrapolation"]
    coeffs = CSV.read(dict["coefficients"], DataFrame)

    scheme = if extrapolation == "none"
        Throw()
    elseif extrapolation == "flat"
        Flat()
    elseif extrapolation == "line"
        Line()
    else
        throw(KeyError("Aerodynamic extrapolation scheme must be 'none', 'flat' or 'line'."))
    end

    if model == "interpolated"
        return from_dict_interpolated(Lref, Sref, XR, coeffs, scheme)
    elseif model == "symmetric"
        return from_dict_symmetric(Lref, Sref, XR, coeffs, scheme)
    else
        throw(KeyError("Aerodynamic model $model not defined."))
    end
end

function from_dict_interpolated(Lref, Sref, XR, coefs, scheme)
    M = coefs.MACH |> unique
    αT = coefs.ALPHA |> unique .|> deg2rad
    ϕA = coefs.PHI |> unique .|> deg2rad
    base = InterpolatedBaseModel(M, αT, ϕA, coefs, scheme)
    deflection = InterpolatedDeflectionModel(M, αT, ϕA, coefs, scheme)
    return ActiveAerodynamics(Lref, Sref, XR, base, deflection)
end

function from_dict_symmetric(Lref, Sref, XR, coefs, scheme)
    M = coefs.MACH |> unique
    αT = coefs.ALPHA |> unique .|> deg2rad
    base = SymmetricBaseModel(M, αT, coefs, scheme)
    deflection = SymmetricDeflectionModel(M, αT, coefs, scheme)
    return ActiveAerodynamics(Lref, Sref, XR, base, deflection)
end

function interp_coeff(M, αT, ϕA, coef, scheme)
    coef_table = permutedims(reshape(coef, (length(ϕA), length(αT), length(M))), (3, 2, 1))
    return extrapolate(interpolate((M, αT, ϕA), coef_table, Gridded(Linear())), scheme)
end

function interp_coeff(M, αT, coef, scheme)
    coef_table = permutedims(reshape(coef, (length(αT), length(M))), (2, 1))
    return extrapolate(interpolate((M, αT), coef_table, Gridded(Linear())), scheme)
end

function from_montecarlo(aed::ActiveAerodynamics, σCA)
    base_mc = deepcopy(aed.base)
    fCA = 1 + randn() * σCA
    base_mc.CAon.itp.coefs .*= fCA
    base_mc.CAoff.itp.coefs .*= fCA
    return ActiveAerodynamics(aed.Lref, aed.Sref, aed.XR, base_mc)
end

end