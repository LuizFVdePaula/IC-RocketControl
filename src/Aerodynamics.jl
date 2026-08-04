module Aerodynamics

export ActiveAerodynamics, aerodynamic_coefficients, getXCP

using CSV
using DataFrames
using Interpolations
using StaticArrays

"""
    struct SymmetricBaseModel

Aerodynamic base coefficients represented as interpolations of (`M`, `αT`) and periodic functions of `ϕA`.

The coefficients are given in aeroballistic coordinate system.
"""
struct SymmetricBaseModel{T}
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

struct SymmetricDeflectionModel{T}
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

struct ActiveAerodynamics
    Lref::Float64
    Sref::Float64
    XR::Float64
    τ::Float64
    base::SymmetricBaseModel
    deflection::SymmetricDeflectionModel
end

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
getXCP(aed::ActiveAerodynamics, M, αT, ΔXCG) = -aed.base.Cm0(M, αT) / aed.base.CN0(M, αT) + ΔXCG

"""
    function aerodynamic_coefficients(aed::ActiveAerodynamics, M, αT, ϕA, ΔXCG, ωBAnorm, δpqr, on::Bool)
    
Obtain aerodynamic coefficents at condition of `M`, `αT`, `ϕA`, `ΔXCG`, `ωBAnorm`, `δpqr`.

Description of inputs:
- `M`: mach number.
- `αT`: total angle of attack.
- `ϕA`: aerodynamic roll angle.
- `ΔXCG`: normalized position difference between reference point and center of gravity (ΔXCG = XR - XCG).
- `ωBAnorm`: normalized angular velocity in body coordinate system.
- `δpqr`: control surface deflections in body coordinate system.
"""
function aerodynamic_coefficients(aed::ActiveAerodynamics, M, αT, ϕA, ΔXCG, ωBAnorm, δpqr, on::Bool)
    sϕA, cϕA = sincos(ϕA)
    TRB = SMatrix{3, 3, Float64, 9}(1, 0, 0, 0, cϕA, sϕA, 0, -sϕA, cϕA)
    ωBAnorm_R = TRB * ωBAnorm
    δpqr_R = TRB * δpqr

    δeff = 0.5 * (abs(δpqr[2]) + abs(δpqr[3]))
    CA = (on ? aed.base.CAon(M, αT) : aed.base.CAoff(M, αT)) + aed.deflection.CAδeff2(M, αT) * δeff^2
    CY = aed.base.CYϕ(M, αT) * sin(4ϕA) + aed.deflection.CYδr(M, αT) * δpqr_R[3]
    CN = aed.base.CN0(M, αT) + aed.base.CNϕ(M, αT) * sin(2ϕA)^2 + aed.deflection.CNδq(M, αT) * δpqr_R[2]
    Cl = aed.base.Clϕ(M, αT) * sin(4ϕA) + aed.base.Clp(M, αT) * ωBAnorm_R[1] + aed.deflection.Clδp(M, αT) * δpqr_R[1]
    Cm = aed.base.Cm0(M, αT) + aed.base.Cmϕ(M, αT) * sin(2ϕA)^2 + aed.base.Cmq(M, αT) * ωBAnorm_R[2] + aed.deflection.Cmδq(M, αT) * δpqr_R[2] - ΔXCG * CN
    Cn = aed.base.Cnϕ(M, αT) * sin(4ϕA) + aed.base.Cnr(M, αT) * ωBAnorm_R[3] + aed.deflection.Cnδr(M, αT) * δpqr_R[3] + ΔXCG * CY

    Fcoef = transpose(TRB) * SVector(CA, CY, CN)
    Mcoef = transpose(TRB) * SVector(Cl, Cm, Cn)
    return (Fcoef, Mcoef)
end

"""
    from_dict(dict::AbstractDict)

Construct `ActiveAerodynamics` from `dict` containing reference length and area, and path to aerodynamic coefficients.

Useful when reading `.json` files.
"""
function from_dict(dict::AbstractDict)
    Lref = dict["reference_length"]
    Sref = dict["reference_area"]
    XR = dict["reference_position"] / Lref
    τ = dict["tau"]
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

    if model == "symmetric"
        return from_dict_symmetric(Lref, Sref, XR, τ, coeffs, scheme)
    else
        throw(KeyError("Aerodynamic model $model not defined."))
    end
end

function from_dict_symmetric(Lref, Sref, XR, τ, coefs, scheme)
    M = coefs.MACH |> unique
    αT = coefs.ALPHA |> unique .|> deg2rad
    base = SymmetricBaseModel(M, αT, coefs, scheme)
    deflection = SymmetricDeflectionModel(M, αT, coefs, scheme)
    return ActiveAerodynamics(Lref, Sref, XR, τ, base, deflection)
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