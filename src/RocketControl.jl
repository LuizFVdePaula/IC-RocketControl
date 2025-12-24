module RocketControl

include("BaseDefs.jl")
include("EnvironmentDefs.jl")
include("StageDefs.jl")
include("Aerodynamics.jl")
include("Propulsion.jl")
include("Structure.jl")
include("Dynamics.jl")
include("RK4Solver.jl")
include("Simulate.jl")

using .BaseDefs, .EnvironmentDefs, .StageDefs, .Simulate

export Environment, environment
export Stage, stage
export simulate

end
