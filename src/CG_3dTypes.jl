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
Point3L = Tuple{Real,Real,Real,Real,Real} # last cord is lift
Sphere3 = Tuple{Real,Real,Real,Real,Real} # x,y,z,q

toPoint3H((x,y,z)::Point3C)::Point3H = (one(x),x,y,z)
lift(x::Float32,y::Float32,z::Float32)::Point3L = (one(x),x,y,z,muladd(widen(z),widen(z), muladd(widen(y),widen(y), widemul(x,x))) ) # promote to double
lift((x,y,z)::Point3C)::Point3L = (one(x),x,y,z, muladd(z,z, muladd(y,y, x*x)) ) # promote to double

dot(u,v) = sum(u.*v)
# cartesian and homogeneous orientation determinants

# cartesian plane equation
function plane(p::Point3C, q::Point3C, r::Point3C)::Plane3
	pq = q.-p
	pr = r.-p
	return (det3(p,pq,pr,1,2,3), -det2(pq,pr,2,3), det2(pq,pr,1,3), -det2(pq,pr,1,2))
end

# computing homogeneous plane equations
plane(p::Point3H,q::Point3H,r::Point3H)::Plane3 = (det3(p,q,r,2,3,4),
													-det3(p,q,r,1,3,4),
													det3(p,q,r,1,2,4),
													-det3(p,q,r,1,2,3))

orientptpl((x,y,z)::Point3C,(pw,px,py,pz)::Plane3) = orientptpl(toPoint3H(x,y,z),(pw,px,py,pz))
orientptpl((w,x,y,z)::Point3H,(pw,px,py,pz)::Plane3) = muladd(px,x,muladd(py,y,muladd(pz,z,pw)))
abovePl(q::Point3,p::Plane3)::Bool = orientptpl(q,p)>zero(p[1])
belowPl(q::Point3,p::Plane3)::Bool = orientptpl(q,p)<zero(p[1])
