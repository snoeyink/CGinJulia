"""
3D Vectors (x, y z), Planes (pw, px, py, pz)
"""

Vector3 = Tuple{Real,Real,Real}
Point3C = Tuple{Real,Real,Real}
Plane3 = Tuple{Real,Real,Real,Real}
Point3H = Tuple{Real,Real,Real,Real}
Point3 = Union{Point3C,Point3H}

toPoint3H((x,y,z)::Point3C)::Point3H = (one(x),x,y,z)

dot((ux,uy,uz)::Tuple{AbstractFloat,AbstractFloat,AbstractFloat},(vx,vy,vz)::Tuple{AbstractFloat,AbstractFloat,AbstractFloat}) = muladd(ux, vx, muladd(uy, vy, (uz * vz))
dot((ux,uy,uz)::Vector3,(vx,vy,vz)::Vector3) = muladd(ux, vx, muladd(uy, vy, (uz * vz))

# cartesian and homogeneous orientation determinants
det2(i::Tuple{Real,Real,Real,Real})::Real = muladd(i[1],i[3],-muladd(i[2],i[4]))
det3(p::Point3C,q::Point3C,r::Point3C)::Real = muladd(p[1],det2(q[2],q[3],r[2],r[3]),muladd(p[2],det2(q[1],q[3],r[1],r[3]),p[3]*det2(q[1],q[2],r[1],r[2])))

# computing homogeneous plane equations
plane(o::Point3C,p::Point3C,q::Point3C,r::Point3C)::Plane3 = (det3(p,q,r),det3(o,q,r),det3(o,p,r),det3(o,p,q))
orientptpl((x,y,z)::Point3C,(pw,px,py,pz)::Plane3) = orientptpl(toPoint3H(x,y,z),(pw,px,py,pz))
orientptpl((w,x,y,z)::Point3H,(pw,px,py,pz)::Plane3) = muladd(px,x,muladd(py,y,muladd(pz,z,pw)))
abovePl(q::Point3,p::Plane3)::Bool = orientptpl(q,p)>zero(p[1])
belowPl(q::Point3,p::Plane3)::Bool = orientptpl(q,p)<zero(p[1])
