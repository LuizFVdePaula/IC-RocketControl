module Simulate

export simulate

using ..BaseDefs: calc_θ, calc_ϕ
using ..RK4Solver
using StaticArrays

function calc_control(sv, t)
    p = sv[11]
    kp = 5
    δp = kp * p
    return [δp, 0.0, 0.0]
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
    N = length(range(trange[begin], trange[end]; step = dt))
    sv_historic = zeros(length(sv₀), N)
    #sv_historic[:, 1] = sv₀

    sv = sv₀
    for (i, t) ∈ enumerate(trange[begin:end-1])
        u = calc_control(sv, t)
        sol = solve(sv, u, stg, env, range(t, t + Ts, step = dt))
        sv_historic[:, 1+(i-1)*n:1+i*n] .= sol
        #u_historic[:, i:i+n] .= u
        sv = sol[:, end]
    end
    return sv_historic
end

end