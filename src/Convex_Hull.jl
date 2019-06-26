function CH(p)
 	#A = []
	D = Dict{Plane3,Bool}()
	F = Tuple{Int64,Int64,Int64}[]
	scatter3d(getindex.(p,2), getindex.(p,3), getindex.(p,4), markershape=:circle,legend=false)
	for i in 1:length(p), j in (i+1):length(p), k in (i+1):length(p)
		(j==k) && continue
		pln = plane(p[i],p[j],p[k])
		if (~(get(D,pln,false)) && ~(all(pln.==0)) && ~(ptsAbove(p,pln)))
			#push!(A,ptsA,ptsB,ptsC)
			#@warn string("""Points: """,ptsA,ptsB,ptsC)
			pts = [p[i];p[j];p[k];p[i]]
			D[pln] = true;
			push!(F,((i-1),(j-1),(k-1)))
			#push!(A,pts)
			#@info string("""Face Values: """,F)
			display(plot3d!(getindex.(pts,2), getindex.(pts,3), getindex.(pts,4), markershape=:circle))
		end
	end
	@warn string("""The length of the vertices array is: """,length(p))
	@info """Please input a filename without the file type (i.e. loremipsum will produce loremipsum.py):"""
	Export("Blender",[string(readline(stdin),""".py"""),p,F])
	#return A
	return p, F
end

function ptsAbove(pts,pl::Plane3)
	for p in pts
		abovePl(p,pl) && return true
	end
		return false
end
