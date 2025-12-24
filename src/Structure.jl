module Structure

using StaticArrays

export StructureSubsystem

"""
    StructureSubsystem

The structures subsystem. Everything in the rocket, except for the propellant, is fixed
and therefore considered part of the structure.

# Fields
- `m`: Rocket mass (without the propellant).
- `Ixx`: Rotational moment of inertia.
- `Iyy`: Longitudinal moment of inertia (without the propellant).
- `xcm`: Center of mass wrt `body system` (without the propellant).
"""
struct StructureSubsystem
    m
    J
    xcm
end

function from_dict(dict::AbstractDict)
    m = dict["mass"]
    Ixx = dict["rotational_moment_of_inertia_x"]
    Iyy = dict["longitudinal_moment_of_inertia_y"]
    Izz = dict["longitudinal_moment_of_inertia_z"]
    Ixy = dict["cross_product_of_inertia_xy"]
    Ixz = dict["cross_product_of_inertia_xz"]
    Iyz = dict["cross_product_of_inertia_yz"]
    xcm = dict["center_of_mass"]
    J = SMatrix{3, 3}([Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz])
    return StructureSubsystem(m, J, xcm)
end

function from_montecarlo(str::StructureSubsystem, σm, σI)
    fm = 1 + randn() * σm
    fI = 1 + randn() * σI
    return StructureSubsystem(str.m * fm, str.J * fI, str.xcm)
end

end