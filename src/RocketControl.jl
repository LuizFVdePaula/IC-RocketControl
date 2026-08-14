module RocketControl

include("BaseDefs.jl")
include("EnvironmentDefs.jl")
include("Aerodynamics.jl")
include("Propulsion.jl")
include("Navigation.jl")
include("Structure.jl")
include("StageDefs.jl")
include("Dynamics.jl")
include("GNC.jl")
include("RK4Solver.jl")
include("Simulate.jl")

using .BaseDefs
using .EnvironmentDefs
using .StageDefs
using .Dynamics
using .GNC
using .Simulate

export Environment, environment, plotinfo
export Stage, stage
export DynamicModel, ControlParameters, continuousmodel, control_fpa, control_pid
export simulate, postprocess, dynamics, ∇, run_monte_carlo

end
