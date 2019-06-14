@testset "Basic 3D Orientation Tests" begin
    w = (0f0,0f0); x = (1,1,2,2); y = (2,2,1,1); z = (1,2,1,2); p = (2f0,2f0,2f0); q = (3,4,3); r = (3f0,1f0,1f0,1f0); s = (2,3,4); t = (-4f0,-1f0,-1f0,-7f0);

    @test plane(x,y,z) == (-3, 3, 3, -3)
    @test dot(q,s) == 30
    @test toPoint3H(q) == (1,q[1],q[2],q[3])
    @test abovePl(t,plane(x,y,z))
    @test !belowPl(toPoint3H(p),plane(x,y,z))

end;
