module CGinJulia
#
export Vector2, Point2, Point2C, Point2H, Line2
export Vector3, Point3, Point3C, Point3H, Plane3, Point3L
export toPoint2H, toPoint3H, lift
export dot, perpdot
export det2, det3, det4
export abovePl, belowPl
export plane
export ccw,leftturn,rightturn
export leftAhead, rightAhead, leftOnAhead, rightOnAhead
export leftofseg, rightofseg, leftward, rightward
export line, orientptln, leftof, rightof

export CH, ptsAbove
export Bl_Export,Export

include("LA_Dets.jl")
include("CG_2dTypes.jl")
include("CG_3dTypes.jl")
include("Export.jl")
include("Blender_Export.jl")
include("Convex_Hull.jl")

end  # module CGinJulia
