module StageDefs

export Stage, stage

import ..Aerodynamics: Aerodynamics, ActiveAerodynamics
import ..Propulsion: Propulsion, SolidPropulsion
import ..Structure: Structure, StructureSubsystem
using JSON

struct Stage
    aed::ActiveAerodynamics
    prp::SolidPropulsion
    str::StructureSubsystem
end

function stage(jsonpath::AbstractString)
    dict = JSON.parse(read(jsonpath, String))
    aed = Aerodynamics.from_dict(dict["Aerodynamics"])
    prp = Propulsion.from_dict(dict["Propulsion"])
    str = Structure.from_dict(dict["STR"])
    return Stage(aed, prp, str)
end

end