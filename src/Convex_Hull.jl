function CH(p)
 	A = []
	F = []
	scatter3d(getindex.(p,2), getindex.(p,3), getindex.(p,4), markershape=:circle,legend=false)
	for ptsA in p[1:length(p)]
		for ptsB in p[1:length(p)]
			for ptsC in p[1:length(p)]
				if (ptsC≠ptsB && ptsC≠ptsA && ptsB≠ptsA)
					pln = plane(ptsA,ptsB,ptsC)
					if (~(all(pln.==0)) && ~(ptsAbove(p,pln)))
						#push!(A,ptsA,ptsB,ptsC)
						@warn string("""Points: """,ptsA,ptsB,ptsC)
						pts = [ptsA;ptsB;ptsC;ptsA]
						push!(F,(findall(f->f==ptsA,p),findall(f->f==ptsB,p),findall(f->f==ptsC,p)))
						push!(A,getindex.(pts,2), getindex.(pts,3), getindex.(pts,4))
						@info string("""Face Values: """,F)
						display(plot3d!(getindex.(pts,2), getindex.(pts,3), getindex.(pts,4), markershape=:circle))
					end
				end
			end
		end
	end
	@warn string("""The length of the vertices array is: """,length(p))
	@info """Please input a filename without the file type (i.e. loremipsum will produce loremipsum.py):"""
	Export("Blender",[string(readline(stdin),""".py"""),p,F])
	return A
end

function ptsAbove(pts,pl::Plane3)
	n = 1
	while n <= length(pts)
		abovePl(pts[n],pl) && return true
		n += 1
	end
		return false
end
