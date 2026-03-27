module StageDefs

export Stage, stage

import ..Aerodynamics: Aerodynamics, ActiveAerodynamics
import ..Propulsion: Propulsion, SolidPropulsion
import ..Structure: Structure, StructureSubsystem
import ..Navigation: Navigation, IMUSensor
using JSON

struct Stage
    aed::ActiveAerodynamics
    prp::SolidPropulsion
    imu::IMUSensor
    str::StructureSubsystem
end

function stage(jsonpath::AbstractString)
    dict = JSON.parse(read(jsonpath, String))
    aed = Aerodynamics.from_dict(dict["Aerodynamics"])
    prp = Propulsion.from_dict(dict["Propulsion"])
    imu = Navigation.from_dict(dict["IMUSensor"])
    str = Structure.from_dict(dict["Structure"])
    return Stage(aed, prp, imu, str)
end

end