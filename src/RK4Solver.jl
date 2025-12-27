module RK4Solver

export solve

import ..Dynamics: dynamics

"""
    nextstatevector(sv, t, dt)

Calculate the next state vector by advancing current `sv` by `dt`.
"""
function nextstatevector(sv, u, stg, env, t, dt)
    k1 = dynamics(sv, u, stg, env, t)
    k2 = dynamics(sv + dt * k1 / 2, u, stg, env, t + dt / 2)
    k3 = dynamics(sv + dt * k2 / 2, u, stg, env, t + dt / 2)
    k4 = dynamics(sv + dt * k3, u, stg, env, t + dt)
    return sv + dt * (k1 + 2k2 + 2k3 + k4) / 6
end

"""
    solve(sv₀, trange::AbstractRange)

Solve (integrate) starting at `sv₀`.
"""
function solve(sv₀, u, stg, env, trange::AbstractRange)
    dt = step(trange)
    hist = zeros(length(sv₀), length(trange))
    sv = sv₀
    for (i, t) ∈ enumerate(trange)
        hist[:, i] = sv
        sv = nextstatevector(sv, u, stg, env, t, dt)
    end
    return hist
end

end