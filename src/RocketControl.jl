module RocketControl

include("BaseDefs.jl")
include("EnvironmentDefs.jl")
include("Aerodynamics.jl")
include("Propulsion.jl")
include("Structure.jl")
include("StageDefs.jl")
include("Dynamics.jl")
include("RK4Solver.jl")
include("Simulate.jl")

using .BaseDefs, .EnvironmentDefs, .StageDefs, .Dynamics, .Simulate

export Environment, environment, plotinfo
export Stage, stage
export simulate, postprocess, dynamics, ∇

end
