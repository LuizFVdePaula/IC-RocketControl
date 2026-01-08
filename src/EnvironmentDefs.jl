module EnvironmentDefs

# TODO
# include ΔT ISA
# wind model
# launch rail: elevation, azimuth and zenith

export gravity, Environment, environment, windspeed, plotinfo

using CairoMakie
using CSV
using DataFrames
using Interpolations
using JSON
using StaticArrays

const gravity = 9.8

"""
    Environment

Contain environmental information, including wind profile.

# Fields
- `wind_north`: Wind noth component as function of altitude.
- `wind_east`: Wind east component as function of altitude.
"""
struct Environment{T}
    wind_north::T
    wind_east::T
end

"""
    Environment(wind_table::DataFrame)

Constructor from wind information provided in `wind_table`.

The `wind_table` must contain a field named `height` (in meters) and wind data expressed in one of the following
formats:
- `vel_north` and `vel_east`: wind velocity components in m/s. A positive `vel_north` indicates wind blowing from south
to north. A positive `vel_east` indicates wind blowing from west to east.
- `speed` and `direction`: wind speed and direction in m/s and degrees, respectively. The direction angle is measured
clockwise from north. A direction of 0° indicates wind blowing from north to south, and a direction of 90° indicates
wind blowing from east to west.
"""
function Environment(wind_table::DataFrame)
    colnames = names(wind_table)
    if "vel_north" ∈ colnames && "vel_east" ∈ colnames
        wind_north = copy(wind_table.vel_north)
        wind_east = copy(wind_table.vel_east)
    elseif "speed" ∈ colnames && "direction" ∈ colnames
        wind_north = -wind_table.speed .* cosd.(wind_table.direction)
        wind_east = -wind_table.speed .* sind.(wind_table.direction)
    else
        KeyError("`wind_table` must contain either fields 'vel_north' and 'vel_east' or 'speed' and 'direction'.")
    end
    return Environment(
        extrapolate(interpolate((wind_table.height,), wind_north, Gridded(Linear())), 0.0),
        extrapolate(interpolate((wind_table.height,), wind_east, Gridded(Linear())), 0.0),
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

function plotinfo(env::Environment)
    hmin, hmax = extrema(env.wind_north.itp.knots[1])
    h = range(hmin, hmax, 1000)
    v_north = env.wind_north(h)
    v_east = env.wind_east(h)
    V = @. sqrt(v_north^2 + v_east^2)

    fig = Figure(size = (900, 600))
    ax_speed = Axis(fig[1, 1], xlabel = "Speed [m/s]", ylabel = "Height AGL [m]")
    lines!(ax_speed, V, h, linewidth = 2)
    ax_direction = PolarAxis(fig[1, 2], title = "Wind direction", direction = -1, theta_0 = -π / 2)
    lines!(ax_direction, atan.(v_east, v_north), h)
    fig
end

end