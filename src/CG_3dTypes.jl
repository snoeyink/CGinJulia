"""
Adapted from 2D Types:
Geometric predicates for 3D Vectors, cartesian points (px,py,pz), homogeneuos points (pw,px,py,pz), and planes
    Julia can inline these.
"""
Vector3 = Tuple{Real,Real,Real}
Point3C = Tuple{Real,Real,Real}
Plane3 = Tuple{Real,Real,Real,Real}
Point3H = Tuple{Real,Real,Real,Real}
Point3 = Union{Point3C,Point3H}

toPoint3H((x,y,z)::Point3C)::Point3H = (one(x),x,y,z)

dot((ux,uy,uz)::Tuple{AbstractFloat,AbstractFloat,AbstractFloat},(vx,vy,vz)::Tuple{AbstractFloat,AbstractFloat,AbstractFloat}) = muladd(ux, vx, muladd(uy, vy, (uz * vz)))
dot((ux,uy,uz)::Vector3,(vx,vy,vz)::Vector3) = muladd(ux, vx, muladd(uy, vy, (uz * vz)))

# cartesian and homogeneous orientation determinants
det2(p,q,a,b)::Real = muladd(p[a],q[b],-p[b]*q[a])
det3(p,q,r,a,b,c)::Real = p[a]*det2(q,r,b,c) -
                          p[b]*det2(q,r,a,c) +
                          p[c]*det2(q,r,a,b)

# computing homogeneous plane equations
plane(o::Point3C,p::Point3C,q::Point3C,r::Point3C)::Plane3 = (det3(p,q,r,2,3,4),
                                                             -det3(o,q,r,1,3,4),
                                                             det3(o,p,r,1,2,4),
                                                             -det3(o,p,q,1,2,3))

orientptpl((x,y,z)::Point3C,(pw,px,py,pz)::Plane3) = orientptpl(toPoint3H(x,y,z),(pw,px,py,pz))
orientptpl((w,x,y,z)::Point3H,(pw,px,py,pz)::Plane3) = muladd(px,x,muladd(py,y,muladd(pz,z,pw)))
abovePl(q::Point3,p::Plane3)::Bool = orientptpl(q,p)>zero(p[1])
belowPl(q::Point3,p::Plane3)::Bool = orientptpl(q,p)<zero(p[1])
