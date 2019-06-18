@testset "Basic 3D Orientation Tests" begin
	w = (0f0,0f0); x = (1,1,2,2); y = (2,2,1,1); z = (1,2,1,2); p = (2f0,2f0,2f0); q = (3,4,3); r = (3f0,1f0,1f0,1f0); s = (2,3,4); t = (-4f0,-1f0,-1f0,-7f0);

	@test plane(x,y,z) == (-3, 3, 3, -3)
	@test dot(q,s) == 30
	@test toPoint3H(q) == (1,q[1],q[2],q[3])
	@test abovePl(t,plane(x,y,z))
	@test !belowPl(toPoint3H(p),plane(x,y,z))

	a = (1f0,3f0,3f0,3f0); b = (1f0,2f0,2f0,2f0); c = (1f0,2f0,2f0,4f0); d = (1f0,4f0,2f0,2f0); e = (1f0,2f0,4f0,2f0); f = (1f0,2f0,4f0,4f0); g = (1f0,4f0,4f0,2f0);
	h = (1f0); i = (1f0);

	pts = [Point3H(1f0, a, b, c) for a in [2f0 4f0], b in [2f0 4f0], c in [2f0 4f0]]

	@test CH([a;b;c;d;e;f;g])

end;
