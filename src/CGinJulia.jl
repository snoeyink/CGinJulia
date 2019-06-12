module CGinJulia
#
export Vector2, Point2, Point2C, Point2H, Line2
export Vector3, Point3, Point3C, Point3H, Line3
export toPoint2H
export dot, perpdot
export ccw,leftturn,rightturn
export leftAhead, rightAhead, leftOnAhead, rightOnAhead
export leftofseg, rightofseg, leftward, rightward
export line, orientptln, leftof, rightof

include("CG_2dTypes.jl")
include("CG_3dTypes.jl")

end  # module CGinJulia
