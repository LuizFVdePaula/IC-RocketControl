module BaseDefs

export calc_ϕ, calc_θ, calc_ψ, calc_q0, calc_q1, calc_q2, calc_q3, rotX, rotY, rotZ, rotXYZ

using StaticArrays

"""
    StateVector{T} <: FieldVector{13, T}

The vehicle current state.

- Position (*x*, *y*, *z*) are given w.r.t. **LVLH** system in NED coordinates.
- Quaternions (*qᵢ*) define **B** orientation w.r.t. **LVLH** system.
- Linear rates (*u*, *v*, *w*) are given w.r.t. **E** frame in **B** system.
- Angular rates (*p*, *q*, *r*) are given w.r.t. **I** frame in **B** system.

# Fields
- `x`: Northing coordinate
- `y`: Easting coordinate
- `z`: Down coordinate
- `qᵢ`: Quaternions (i ∈ {0, 1, 2, 3})
- `u`: Forward speed
- `v`: Side speed
- `w`: Normal speed
- `p`: Roll rate
- `q`: Pitch rate
- `r`: Yaw rate
"""
mutable struct StateVector{T} <: FieldVector{13, T}
    φ::T
    λ::T
    h::T
    q0::T
    q1::T
    q2::T
    q3::T
    u::T
    v::T
    w::T
    p::T
    q::T
    r::T
end

calc_ϕ(q0, q1, q2, q3) = atan(2 * (q2 * q3 + q0 * q1), q0^2 - q1^2 - q2^2 - q3^2)
calc_θ(q0, q1, q2, q3) = asin(-2 * (q1 * q3 - q0 * q2))
calc_ψ(q0, q1, q2, q3) = atan(2 * (q1 * q2 + q0 * q3), q0^2 + q1^2 - q2^2 - q3^2)

calc_ϕ(t13, t23, t33) = acos(t33 / cos(asin(-t13))) * sign(t23)
calc_θ(t13) = asin(-t13)
calc_ψ(t11, t12, t13) = acos(t11 / cos(asin(-t13))) * sign(t12)

calc_ϕ(sv::StateVector) = calc_ϕ(sv.q0, sv.q1, sv.q2, sv.q3)
calc_θ(sv::StateVector) = calc_θ(sv.q0, sv.q1, sv.q2, sv.q3)
calc_ψ(sv::StateVector) = calc_ψ(sv.q0, sv.q1, sv.q2, sv.q3)

calc_q0(ϕ, θ, ψ) = cos(ψ / 2) * cos(θ / 2) * cos(ϕ / 2) + sin(ψ / 2) * sin(θ / 2) * sin(ϕ / 2)
calc_q1(ϕ, θ, ψ) = cos(ψ / 2) * cos(θ / 2) * sin(ϕ / 2) - sin(ψ / 2) * sin(θ / 2) * cos(ϕ / 2)
calc_q2(ϕ, θ, ψ) = cos(ψ / 2) * sin(θ / 2) * cos(ϕ / 2) + sin(ψ / 2) * cos(θ / 2) * sin(ϕ / 2)
calc_q3(ϕ, θ, ψ) = sin(ψ / 2) * cos(θ / 2) * cos(ϕ / 2) - cos(ψ / 2) * sin(θ / 2) * sin(ϕ / 2)

rotX(ϕ) = SMatrix{3, 3}([1 0 0; 0 cos(ϕ) sin(ϕ); 0 -sin(ϕ) cos(ϕ)])
rotY(θ) = SMatrix{3, 3}([cos(θ) 0 -sin(θ); 0 1 0; sin(θ) 0 cos(θ)])
rotZ(ψ) = SMatrix{3, 3}([cos(ψ) sin(ψ) 0; -sin(ψ) cos(ψ) 0; 0 0 1])

function rotXYZ(q0, q1, q2, q3)
    sx, sy, sz = 2 * q0 * q1,   2 * q0 * q2,    2 * q0 * q3
    xx, xy, xz = 2 * q1^2,      2 * q1 * q2,    2 * q1 * q3
    yy, yz, zz = 2 * q2^2,      2 * q2 * q3,    2 * q3^2
    r = SMatrix{3, 3}([
        1 - (yy + zz)   xy + sz         xz - sy
        xy - sz         1 - (xx + zz)   yz + sx
        xz + sy         yz - sx         1 - (xx + yy)
    ])
    return r
end

rotXYZ(sv::StateVector) = rotXYZ(sv.q0, sv.q1, sv.q2, sv.q3)

end