@testset "Basic 3D Orientation Tests" begin
    w = (0f0,0f0); x = (1,1,2,2); y = (2,2,1,1); z = (1,2,1,2); p = (2,2,2); q = (3,4,3); r = (3f0,1f0,1f0,1f0); s = (2,3,4);

    @test plane(r,x,y,z)

end;
