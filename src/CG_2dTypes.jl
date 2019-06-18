""" geometric predicates for 2vectors, cartesian points (px,py), homogeneuos points (pw,px,py), lines (lw,lx,ly), and segments x[1:2],y[1:2].
	Julia can inline these.
"""
Vector2 = Tuple{Real,Real}
Point2C = Tuple{Real,Real}
Line2 = Tuple{Real,Real,Real}
Point2H = Tuple{Real,Real,Real}
Point2 = Union{Point2C,Point2H}

toPoint2H((x,y)::Point2C)::Point2H = (one(x),x,y)

# some special handling for floats
dot((ux,uy)::Tuple{AbstractFloat,AbstractFloat}, (vx,vy)::Tuple{AbstractFloat,AbstractFloat}) = muladd(ux, vx, 1.0uy*vy)
dot((ux,uy)::Vector2, (vx,vy)::Vector2) = muladd(ux, vx, uy*vy) # dot product for 2d vectors. (consider muladd)
perpdot((ux,uy)::Tuple{AbstractFloat,AbstractFloat}, (vx,vy)::Tuple{AbstractFloat,AbstractFloat}) = muladd(ux, vy, -1.0uy*vx)
perpdot((ux,uy)::Vector2, (vx,vy)::Vector2) = muladd(ux, vy, -uy*vx) # dot product of ccw perp(u) with v. (consider muladd)

# basic cartesian and homogeneous orientation determinants
ccw(p::Point2C,q::Point2C,r::Point2C) = perpdot(q.-p,r.-p) # orientation determinant
ccw((pw,px,py)::Point2H,(qw,qx,qy)::Point2H,(rw,rx,ry)::Point2H) = pw*perpdot((qx,qy),(rx,ry)) -qw*perpdot((px,py),(rx,ry)) + rw*perpdot((px,py),(qx,qy))# orientation determinant
ccw(p::Point2C,q::Point2C,(rw,rx,ry)::Point2H) = perpdot(q.-p,(rx,ry))+rw*perpdot(p,q.-p) # orientation determinant

# mixed cartesian and homogeneuos
ccw(p::Point2C,q::Point2H,r::Point2H) = ccw(toPoint2H(p),q,r)
ccw(p::Point2H,q::Point2C,r::Point2H) = ccw(p,toPoint2H(q),r)
ccw(p::Point2H,q::Point2H,r::Point2C) = ccw(p,q,toPoint2H(r))
ccw(p::Point2H,q::Point2C,r::Point2C) = ccw(q,r,p)
ccw(p::Point2C,q::Point2H,r::Point2C) = ccw(r,p,q)

## Geometric predicates based on ccw
leftturn(p::Point2,q::Point2,r::Point2)::Bool = ccw(p,q,r)>zero(p[1]) # p to q to r is a left turn
rightturn(p::Point2,q::Point2,r::Point2)::Bool = ccw(p,q,r)<zero(p[1]) # p to q to r is a right turn

## leftAhead halps solve on-the-llne cases consistent with infinitesimal rotation
leftAhead(p::Point2C,q::Point2C,r::Point2C)::Bool = (d = ccw(p,q,r)) == zero(d) ? dot(p.-q,r.-q)<zero(d) : d>zero(d) # should assert ≉ 0 for precision checking
rightAhead(p::Point2C,q::Point2C,r::Point2C)::Bool = (d = ccw(p,q,r)) == zero(d) ? dot(p.-q,r.-q)<zero(d) : d<zero(d) # should assert ≉ 0 for precision checking
leftOnAhead(p::Point2C,q::Point2C,r::Point2C)::Bool = (d = ccw(p,q,r)) == zero(d) ? dot(q.-p,r.-p)>zero(d) : d>zero(d) # should assert ≉ 0 for precision checking
rightOnAhead(p::Point2C,q::Point2C,r::Point2C)::Bool = (d = ccw(p,q,r)) == zero(d) ? dot(q.-p,r.-p)>zero(d) : d<zero(d) # should assert ≉ 0 for precision checking

## orientation with drected segments
#  This should probably be expanded. 
leftofseg(q::Point2, xs, ys)::Bool = @inbounds leftturn((xs[1],ys[1]),(xs[2],ys[2]), q) # point q left of segment defined by pairs of xs and ys
rightofseg(q::Point2, xs, ys)::Bool = @inbounds rightturn((xs[1],ys[1]),(xs[2],ys[2]), q) # point q right of segment defined by pairs of xs and ys

leftward(s0::Point2C,s1::Point2C, t0::Point2C,t1::Point2C)::Bool = perpdot(t1.-t0, s1.-s0)>zero(s0[1])   # vector s0->s1 is in the direction from right to left of the line t0->t1.
rightward(s0::Point2C,s1::Point2C, t0::Point2C,t1::Point2C)::Bool = perpdot(t1.-t0, s1.-s0)<zero(s0[1])   # vector s0->s1 is in the direction from left to right of the line t0->t1.

# computing homogeneous line equations
line(p::Point2C,q::Point2C)::Line2 = (perpdot(p,q), last(p)-last(q), first(q)-first(p)) # line equation
orientptln((x,y)::Point2C, (lw,lx,ly)::Line2) = muladd(y,ly,muladd(x,lx,lw))
orientptln((w,x,y)::Point2H, (lw,lx,ly)::Line2) = muladd(y,ly,muladd(x,lx,w*lw))
leftof(q::Point2, l::Line2)::Bool = orientptln(q,l)>zero(l[1]) # point left of (consistent with) oriented line
rightof(q::Point2, l::Line2)::Bool = orientptln(q,l)<zero(l[1])
