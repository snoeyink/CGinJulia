using Plots

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
plane(p::Point3H,q::Point3H,r::Point3H)::Plane3 = (det3(p,q,r,2,3,4),
                                                  -det3(p,q,r,1,3,4),
                                                  det3(p,q,r,1,2,4),
                                                  -det3(p,q,r,1,2,3))

orientptpl((x,y,z)::Point3C,(pw,px,py,pz)::Plane3) = orientptpl(toPoint3H(x,y,z),(pw,px,py,pz))
orientptpl((w,x,y,z)::Point3H,(pw,px,py,pz)::Plane3) = muladd(px,x,muladd(py,y,muladd(pz,z,pw)))
abovePl(q::Point3,p::Plane3)::Bool = orientptpl(q,p)>zero(p[1])
belowPl(q::Point3,p::Plane3)::Bool = orientptpl(q,p)<zero(p[1])

function CH(p)
  A = []
  F = []
  scatter3d(getindex.(p,2), getindex.(p,3), getindex.(p,4), markershape=:circle,legend=false)
  for ptsA in p
    for ptsB in p[1:length(p)]
      for ptsC in p[1:length(p)]
        (ptsC==ptsB || ptsC==ptsA || ptsB==ptsA) && continue
          pln = plane(ptsA,ptsB,ptsC)
          if (~(all(pln.==0)) && ~(ptsAbove(p,pln)))
            #push!(A,ptsA,ptsB,ptsC)
            pts = [ptsA;ptsB;ptsC;ptsA]
            push!(A,getindex.(pts,2), getindex.(pts,3), getindex.(pts,4))
            display(plot3d!(getindex.(pts,2), getindex.(pts,3), getindex.(pts,4), markershape=:circle))
        end
      end
    end
  end
  Bl_Export("test1.py",A,F)
  return A
end

function ptsAbove(pts,pl::Plane3)
  n = 1
    while n <= length(pts)
      abovePl(pts[n],pl) && return true
      n += 1
    end
      return false
end
