module FloatTiming

using StaticArrays, DoubleFloats, Test, BenchmarkTools

widemult(a::Float64, b::Real) = Double64(a)*Double64(b)
widemult(a::Union{Float32,Int}, b::Float64) = Double64(a)*Double64(b)
widemult(a, b) = widemul(a,b)

""" det2(p,q, a,b)
2x2 same precision
"""
det2(p,q, a,b) = (p[a]*q[b]) - (p[b]*q[a])
det2f(p,q, a,b) = muladd(p[a],q[b], -(p[b]*q[a]))
det2w(p,q, a,b) = widemult(p[a],q[b]) - widemult(p[b],q[a]) # 2 mult, add, neg
det2wm(p,q, a,b) = muladd(p[a],q[b], -widemult(p[b],q[a])) # 2 mult, add, neg
@inline det2i(p,q, a,b) =  @inbounds (p[a]*q[b] - p[b]*q[a])
@inline det2fi(p,q, a,b) = @inbounds muladd(p[a],q[b], -(p[b]*q[a]))
@inline det2wi(p,q, a,b) = @inbounds widemult(p[a],q[b]) - widemult(p[b],q[a]) # 2 mult, add, neg
@inline det2mi(p,q, a,b) = @inbounds muladd(p[a],q[b], -widemult(p[b],q[a])) # 2 mult, add, neg

function test(alg,pts)
	y = @MVector [alg(pts[1],pts[2],1,2), alg(pts[1],pts[2],1,3), alg(pts[1],pts[2],2,3)]
	for p in pts, q in pts
		y[1] += alg(p,q,1,2)
		y[2] += alg(p,q,1,3)
		y[3] += alg(p,q,2,3)
	end
	y
end

det2algs = [det2, det2f, det2w, det2wm, det2i,det2fi,det2wi,det2mi]
pts = [rand(Float32,3) for i in 1:20]

for alg in det2algs
	@time test(alg, pts)
end

det3(p,q,r, a,b,c) = muladd(p[a], det2(q,r,b,c), # 9 mult, 5 add, 3 neg,
					 muladd(p[b], det2(r,q,a,c),
						    p[c]*det2(q,r,a,b)))
det4ck(p,q,r,s, a,b,c,d) = muladd(p[a], det3(q,r,s, b,c,d), # 40 mult, 23 add, 12 neg
						 muladd(p[b], det3(q,r,s, a,d,c),
						 muladd(p[c], det3(q,r,s, a,b,d),
								p[d]*det3(q,r,s, a,c,b))))

# 2x2 minors is better
det4(p,q,r,s, a,b,c,d) = Double64(det2(p,q, a,b))*Double64(det2(r,s, c,d)) - # 30 mult, 17 add, 12 neg
						 Double64(det2(p,q, a,c))*Double64(det2(r,s, b,d)) +
						 Double64(det2(p,q, a,d))*Double64(det2(r,s, b,c)) +
						 Double64(det2(p,q, b,c))*Double64(det2(r,s, a,d)) -
						 Double64(det2(p,q, b,d))*Double64(det2(r,s, a,c)) +
						 Double64(det2(p,q, c,d))*Double64(det2(r,s, a,b))

Plucker3L(p,q) = (q.-p, det2(p,q, 1,2), det2(p,q, 1,3), det2(p,q, 1,4), det2(p,q, 2,3), det2(p,q, 2,4), det(p,q, 3,4))

det4dbl(p,q,r,s, a,b,c,d) = muladd(det2(p,q, a,b), det2(r,s, c,d), # 30 mult, 17 add, 12 neg
						 muladd(det2(p,q, c,a), det2(r,s, b,d),
						 muladd(det2(p,q, a,d), det2(r,s, b,c),
						 muladd(det2(p,q, b,c), det2(r,s, a,d),
						 muladd(det2(p,q, d,b), det2(r,s, a,c), det2(p,q, c,d))*det2(r,s, a,b)) ) ) )





end
