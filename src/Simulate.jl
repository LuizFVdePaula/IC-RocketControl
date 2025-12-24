module Simulate

export simulate

using ..RK4Solver
using StaticArrays

function calc_control(sv)
    return zeros(3)
end

"""
1 -> sv₀ (t = 0.0)
solve resolve de [1:11] (t = 0.0 até 1.0)
append sol[2:11] (t = 0.1 até 1.0)
sv = sol[end]

"""
function simulate(stg, env, sv₀, trange)
    n = 10
    Ts = step(trange)
    dt = Ts / n
    N = length(range(tstart, tstop; step = dt))
    sv_historic = zeros(length(sv₀), N)
    sv_historic[:, 1] = sv₀

    sv = sv₀
    for (i, t) ∈ enumerate(trange)
        u = calc_control(sv)
        sol = solve(sv, u, stg, env, range(t, t + Ts, step = dt))
        sv_historic[:, i:i+n] = sol[2:end]
        #u_historic[:, i:i+n] .= u
        sv = sol[end]
    end
    return sv_historic
end

end