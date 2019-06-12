@testset "basic 2d orientation tests" begin
    o = (0f0,0f0); x = (1,0); y = (0,1); p = (2,2); q = (3,4); r = (1f0,1f0);
    @test ccw(o,x,p) == ccw(o,x,toPointWXY(p))
    @test ccw(o,x,p) == ccw(o,toPointWXY(x),(p))
    @test ccw(o,x,p) == ccw(toPointWXY(o),x,(p))
    @test ccw(o,x,p) == ccw(x,p,toPointWXY(o))
    @test ccw(o,x,p) == ccw(toPointWXY(o),toPointWXY(x),toPointWXY(p))
    @test ccw(o,x,p) == ccw(toPointWXY(o),toPointWXY(x),(p))
    @test ccw(o,x,p) == ccw(toPointWXY(o),(x),toPointWXY(p))
    @test ccw(o,x,p) == ccw((o),toPointWXY(x),toPointWXY(p))
    @test 0 == ccw(o,r,toPointWXY(p))
    @test 0 == ccw(r,p,toPointWXY(o))
    @test 0 == ccw((o),toPointWXY(r),toPointWXY(p))
    @test 0 == ccw(toPointWXY(o),toPointWXY(r),(p))
    @test 0 == ccw(toPointWXY(o),toPointWXY(r),toPointWXY(p))

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
    @test rightof(toPointWXY(p),line(r,q))
    @test leftof(toPointWXY(r),line(p,q))
    @test !leftof(toPointWXY(o), line(p,r))
    @test !rightof(toPointWXY(o), line(p,r))
end;
