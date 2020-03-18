struct Edge
	from::Int
	to::Int
	cap
end

function Find!(n::Int, dsu::Array{Int64,1})
	if (dsu[n]==n)
		return n
	else
		k=dsu[n];
		z=Find!(k,dsu);
		dsu[n]=z;
		return z
	end
end

function Union!(x::Int,y::Int,dsu::Array{Int64,1},sz::Array{Int64,1})
	X=Find!(x,dsu);
	Y=Find!(y,dsu);
	if (X==Y)
		return nothing;
	else
		if (sz[X]<sz[Y])
			sz[Y]+=sz[X];
			dsu[X]=Y;
		else
			sz[X]+=sz[Y];
			dsu[Y]=X;
		end
	end
	return nothing;
end


function get_bottleneck!( s::Int,t::Int,n::Int , m::Int , edge_list ::Array{Edge,1} )
	sort!(edge_list,by = x-> -x.cap);
	sz=fill(1,n);
	dsu=[iter for iter in 1:n];

	for i in edge_list
		Union!(i.from,i.to,dsu,sz);
		if (dsu[s]==dsu[t])
			return i.cap;
		end
	end
	return 0;
end 

struct Unconnected_Graph <: Exception end

function change_edges!(s::Int,t::Int,n::Int,m::Int,edge_list::Array{Edge,1},eps)
	B=get_bottleneck!(s,t,n,m,edge_list);
	if (B==0)
		throw(Unconnected_Graph);
	end
	lb=eps*B/(2*m);
	ub=m*B;
	cnt=m;
	for i in 1:m
		if (edge_list[i].cap>ub)
			edge_list[i]=Edge(edge_list[i].from,edge_list[i].to,ub);
		end
		if (edge_list[i].cap<lb)
			cnt=i-1;
			break;
		end
	end
	edge_list[1:cnt]
end
