@testset "Basic 2D Orientation Tests" begin
    o = (0f0,0f0); x = (1,0); y = (0,1); p = (2,2); q = (3,4); r = (1f0,1f0);
    @test ccw(o,x,p) == ccw(o,x,toPoint2H(p))
    @test ccw(o,x,p) == ccw(o,toPoint2H(x),(p))
    @test ccw(o,x,p) == ccw(toPoint2H(o),x,(p))
    @test ccw(o,x,p) == ccw(x,p,toPoint2H(o))
    @test ccw(o,x,p) == ccw(toPoint2H(o),toPoint2H(x),toPoint2H(p))
    @test ccw(o,x,p) == ccw(toPoint2H(o),toPoint2H(x),(p))
    @test ccw(o,x,p) == ccw(toPoint2H(o),(x),toPoint2H(p))
    @test ccw(o,x,p) == ccw((o),toPoint2H(x),toPoint2H(p))
    @test 0 == ccw(o,r,toPoint2H(p))
    @test 0 == ccw(r,p,toPoint2H(o))
    @test 0 == ccw((o),toPoint2H(r),toPoint2H(p))
    @test 0 == ccw(toPoint2H(o),toPoint2H(r),(p))
    @test 0 == ccw(toPoint2H(o),toPoint2H(r),toPoint2H(p))

    @test leftturn(o,x,p)
    @test !rightAhead(o,x,p)
    @test rightturn(o,y,q)
    @test !leftAhead(o,y,q)
    @test !leftturn(x,p,q)
    @test leftAhead(x,p,q) # colinear
    @test rightAhead(x,p,q) # colinear
    @test !leftAhead(o,p,r) # colinear
    @test leftOnAhead(y,o,q)
    @test leftOnAhead(o,p,r) # colinear
    @test rightOnAhead(o,p,r) # colinear
    @test rightOnAhead(x,p,q) # colinear
    @test !rightOnAhead(o,x,p) # colinear
    @test leftofseg(p,y,x)  # query point with segment
    @test rightofseg(p,x,y)
    @test leftward(p,r, o,y)
    @test !rightward(p,r, o,y)
    @test !leftward(p,r, o,r) # parallel
    @test !rightward(p,r, o,r) # parallel
    @test line(o,y) == (0,-1,0)
    @test line(y,r) == (-1,0,1)
    @test leftturn(line(o,x), line(x,p), line(p,o))
    @test rightof(p,line(r,q))
    @test leftof(r,line(p,q))
    @test !leftof(o, line(p,r))
    @test !rightof(o, line(p,r))
    @test rightof(toPoint2H(p),line(r,q))
    @test leftof(toPoint2H(r),line(p,q))
    @test !leftof(toPoint2H(o), line(p,r))
    @test !rightof(toPoint2H(o), line(p,r))
end;
