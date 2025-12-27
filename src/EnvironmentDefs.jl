module EnvironmentDefs

export Environment, environment, windspeed

using CSV
using DataFrames
using Interpolations
using JSON
using StaticArrays

"""
    Environment

Contain environmental information, including launch rail and wind profile.

# Fields
- `wind_north`: Wind noth component as function of altitude.
- `wind_east`: Wind east component as function of altitude.
- `g`: Gravitational acceleration.
"""
struct Environment
    wind_north
    wind_east
    g
end

"""
    Environment(wind_table::DataFrame)

Constructor from `wind_table`, which must have `height`, `vel_north` and `vel_east` fields.
"""
function Environment(wind_table::DataFrame)
    return Environment(
        linear_interpolation(wind_table.height, wind_table.vel_north, extrapolation_bc = 0),
        linear_interpolation(wind_table.height, wind_table.vel_east, extrapolation_bc = 0),
        9.8,
    )
end

"""
    environment(jsonpath::AbstractString)

Create a `Environment` from a **.json** file containing environment data.
"""
function environment(jsonpath::AbstractString)
    dict = JSON.parse(read(jsonpath, String))
    wind_table = CSV.read(dict["wind"], DataFrame)
    return Environment(wind_table)
end

windspeed(env::Environment, h) = SVector(env.wind_north(h), env.wind_east(h), 0.0)

end