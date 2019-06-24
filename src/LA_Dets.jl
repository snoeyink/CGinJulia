# cartesian and homogeneous orientation determinants
det2(p,q,a,b)::Real = muladd(p[a],q[b],-p[b]*q[a])
det3(p,q,r,a,b,c)::Real = p[a]*det2(q,r,b,c) -
							p[b]*det2(q,r,a,c) +
							p[c]*det2(q,r,a,b)
