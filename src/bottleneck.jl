mutable struct Edge
	from		# src vertex id {1..n}
	to			# dest vertex id {1..n}
	cap			# capacity
	wt			# weight (used in mult wt)
	cond		# conductance (used in oracle)
	flow		# current flow (filled by oracle)
	acc_flow	# accumulated flow (used in mult wts)
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
	# disjoint set union with path compression and size heuristic O(alpha(n)) per query amortized
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


# In all functions s:source, t:sink, n:number of vertices, m:number of edges, edge_list: array of edges

function get_bottleneck!( s::Int,t::Int,n::Int , m::Int , edge_list ::Array{Edge,1} )
	# returns bottleneck b and also sorts edge_list (pass by reference) in descending order by capacity
	# complexity O(mlogm) dominated by sorting
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
	# Does changes in the edge_list also!! Returns final edge list
	B=get_bottleneck!(s,t,n,m,edge_list);
	if (B==0)
		throw(Unconnected_Graph);
	end
	lb=eps*B/(2*m);
	ub=m*B;
	cnt=m;
	for i in 1:m
		if (edge_list[i].cap>ub)
			edge_list[i].cap=ub;
		end
		if (edge_list[i].cap<lb)
			cnt=i-1;
			break;
		end
	end
	(B,edge_list[1:cnt])
end

