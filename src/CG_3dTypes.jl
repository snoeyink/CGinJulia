"""

"""

Vector3 = Tuple{Real,Real,Real}
Point3C = Tuple{Real,Real,Real}
Line3 = Tuple{Real,Real,Real,Real}
Point3H = Tuple{Real,Real,Real,Real}
Point3 = Union{Point3C,Point3H}

toPoint3H((x,y,z)::Point3C)::Point3H = (one(x),x,y,z)

dot((ux,uy,uz)::Tuple{AbstractFloat,AbstractFloat,AbstractFloat},(vx,vy,vz)::Tuple{AbstractFloat,AbstractFloat,AbstractFloat}) = (ux * vx) + (uy * vy) + (uz * vz)
dot((ux,uy,uz)::Vector3,(vx,vy,vz)::Vector3) = (ux * vx) + (uy * vy) + (uz * vz)
"""
Above: not getting correct answer, probably a flaw in logic/understanding. Will continue.
"""
