module CGinJulia
#
export Vector2, Point2, Point2C, Point2H, Line2
export Vector3, Point3, Point3C, Point3H, Plane3
export toPoint2H, toPoint3H
export dot, perpdot
export det2, det3
export abovePl, belowPl
export plane
export CH, ptsAbove
export ccw,leftturn,rightturn
export leftAhead, rightAhead, leftOnAhead, rightOnAhead
export leftofseg, rightofseg, leftward, rightward
export line, orientptln, leftof, rightof
export Bl_Export,Export

include("CG_2dTypes.jl")
include("CG_3dTypes.jl")
include("Export.jl")
include("Blender_Export.jl")
include("Convex_Hull.jl")
include("LA_Dets.jl")

end  # module CGinJulia
