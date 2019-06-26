# cartesian and homogeneous orientation determinants
using DoubleFloats
#  iszero, isnonzero, isone                 #  value == 0, value != 0, value == 1
#  ispositive, isnegative,                  #  value >  0, value <  0

# I don't want Float64 to widen to BigFloat, but to Double64
widemult(a::Real, b::Float64) = Double64(a)*Double64(b)
widemult(a::Float64, b::Real) = Double64(a)*Double64(b)
widemult(a, b) = widemul(a,b)

det2(p,q, a,b) = muladd(p[a],q[b], widemult(-p[b],q[a]))
det3(p,q,r, a,b,c) = muladd(p[a],det2(q,r,b,c),
						muladd(p[b], det2(r,q,a,c),
							   p[c]*det2(q,r,a,b)))
det4(p,q,r,s, a,b,c,d) = muladd(p[a],det3(q,r,s, b,c,d),
							muladd(p[b], det3(r,q,s, a,d,c),
							muladd(p[c], det3(r,q,s, a,b,d),
								   p[d]*det3(q,r,s, a,c,b))))
