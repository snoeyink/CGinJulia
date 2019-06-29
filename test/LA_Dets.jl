import CGinJulia.widemult, CGinJulia.det2, CGinJulia.det3, CGinJulia.det4, CGinJulia.det4ck

## I don't want Float64 to widen to BigFloat, but to Double64
@testset "wide multiplication and DoubleFloats" begin
	f = 1f0; g = 2f0; h = 2f0^22-1
	df = 1.0; dg = 2.0; dh = 2.0^52-1
	@test(typeof(h) == Float32)
	@test(typeof(widemult(f,g))==Float64)
	@test(typeof(widemult(df,dg))==Double64)
	@test(widemult(g,h)+2 == 2.0^23)
	@test(widemult(dg,dh)+2 == 2.0^53)
end
##

@testset "determinants" begin
	e1 = (1,0,0,0); e2 = (0,1,0,0); e3 = (0,0,1,0); e4 = (0,0,0,1);
	r = (1,2,3,5); s = (7,11,13,17); t = (1, 2, 4, 8); u = (1, 10, 100, 1000)
	@test(det2(e1,e2, 1,2) == 1)
	@test(det3(e1,e2,e3, 1,2,3) == 1)
	@test(det4(e1,e2,e3,e4, 1,2,3,4) == 1)
	@test(det2(r,s, 1,2) == -3)
	@test(det3(r,s,t, 1,2,3) == -3)
	@test(det4(r,s,t,u, 1,2,3,4) == -2160)
	@test(det4(r,s,t,u, 1,2,3,4) == det4ck(r,s,t,u, 1,2,3,4))
end
