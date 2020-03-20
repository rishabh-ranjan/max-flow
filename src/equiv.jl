include("struct_defs.jl")

using DataStructures

# returns maximum bottleneck and
# sorts edges in descending order by capacity
# complexity O(mlogm) dominated by sorting
function get_max_bottleneck!(g::Graph)
	sort!(g.edges, by = x -> -x.cap)
	dsu = IntDisjointSets(g.n)
	for e in g.edges
		IntDisjointSets.union!(dsu, e.from, e.to)
		if IntDisjointSets.in_same_sets(dsu, g.s, g.t)
			return e.cap
		end
	end
	# here iff s and t are disconnected in graph
	return 0.0
end 

# modifies graph to convert into an equivalent problem
# but with polynomially bounded ratio of edge capacities
function equiv_graph!(g::Graph)
	# TODO: customize eps/2.0
	b = get_max_bottleneck!(g)
	lb = g.eps*b/(2*g.m)
	ub = g.m*b
	cnt = g.m
	for e in edges
		if e.cap > ub
			e.cap = ub
		elseif e.cap < lb
			cnt = i - 1
			break
		end
	end
	g.eps /= 2.0
	g.m = cnt
	g.edges = g.edges[1:cnt]
end

