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

The solution contains every state vector ranging `trange`.
"""
function solve(sv₀, u, stg, env, trange::AbstractRange)
    dt = step(trange)
    N = length(trange)
    hist = Matrix{eltype(sv₀)}(undef, length(sv₀), N)
    sv = sv₀
    for i ∈ 1:(N-1)
        hist[:, i] = sv
        t = trange[i]
        sv = nextstatevector(sv, u, stg, env, t, dt)
    end
    hist[:, N] = sv
    return hist
end

end