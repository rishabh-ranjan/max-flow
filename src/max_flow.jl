# types have been specified in the hopes of compiler optimization

# internal representation of edges
mutable struct Edge
	from::Int64			# src vertex id {1..n}
	to::Int64			# dest vertex id {1..n}
	cap::Float64		# capacity
	wt::Float64			# weight (used in mult wt)
	cond::Float64		# conductance (used in oracle)
	flow::Float64		# current flow (filled by oracle)
	acc_flow::Float64	# accumulated flow (used in mult wts)
end

# collection of parameters describing the approx max-flow problem
# note that we solve for undirected graphs
# even though the representation is directed
# `eps' and `edges' may be modified
mutable struct Graph
	eps::Float64		# epsilon for approximation
	n::Int64			# number of vertices
	m::Int64			# number of edges
	s::Int64			# source vertex
	t::Int64			# sink vertex
	edges::Array{Edge, 1}	# edges in the graph
	b::Float64			# max bottleneck in the graph
end

using DataStructures

# sets maximum bottleneck of the graph and
# sorts edges in descending order by capacity
# return type void
# complexity O(mlogm) dominated by sorting
function set_max_bottleneck!(g::Graph)
	sort!(g.edges, by = x -> -x.cap)
	dsu = IntDisjointSets(g.n)
	for e in g.edges
		IntDisjointSets.union!(dsu, e.from, e.to)
		if IntDisjointSets.in_same_sets(dsu, g.s, g.t)
			g.b = e.cap
			return
		end
	end
	# here iff s and t are disconnected in graph
	g.b = 0
end 

# modifies graph to convert into an equivalent problem
# but with polynomially bounded ratio of edge capacities
function equiv_graph!(g::Graph)
	# TODO: customize eps/2.0
	lb = g.eps*g.b/(2*g.m)
	ub = g.m*g.b
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

# (eps, 3*sqrt(m/eps)) oracle
# modifies the flows of the edges in-place
# returns false for "fail", else returns true
# Note: flows get corrupted (these are no more needed anyway)
# even when false is returned
function oracle!(g::Graph, flow_val::Int64)

	# sum of wts
	ws = 0.0
	for e in g.edges
		ws += e.wt
	end
	# calculation of conductances
	for e in g.edges
		e.cond = e.cap^2 / (e.wt + (g.eps*ws)/(3*g.m))
	end

	# construction of laplacian matrix
	# sparse matrix has to be used as
	# simple matrix construction and use would be O(n^2)
	# dominating the intended complexity
	lap = SparseArrays.spzeros(g.n, g.n)
	for e in g.edges
		lap[e.from, e.to] = -e.cond
		lap[e.to, e.from] = -e.cond
		lap[e.from, e.from] += e.cond
		lap[e.to, e.to] += e.cond
	end

	# construction of chi_s,t
	# TODO: this need not be recomputed, move it out
	chi = SparseArrays.spzeros(g.n)
	chi[g.s] = 1
	chi[g.t] = -1

	# TODO: pass the del = eps/3 parameter
	solver = Laplacians.KMPSDDMSolver(lap)
	# vertex potentials
	phi = solver(chi)

	for e in g.edges
		e.flow = e.cond * (phi[e.to] - phi[e.from]) * flow_val
	end

	# calculate energy of the flow
	eng = 0.0
	for e in g.edges
		eng += e.flow^2 / e.cond
	end

	# TODO: move it inside the loop for micro-optimization
	return eng <= (1+g.eps)*ws
end

# uses the multiplicative weights method to
# either compute a flow with value within O(g.eps) approximation of flow_val or "fail"
# the flow is modified in the graph itself
# false is returned on failure and true on success
# the flows may get corrupted even when false is returned
function mult_wt!(g::Graph, flow_val::Float64, oracle!::Function, rho::Float64)
	for e in g.edges
		e.wt = 1.0
		e.acc_flow = 0.0
	end
	# number of iterations
	n = ceil(Int, 2 * rho * max(log(m),1) / eps^2)
	for i in 1:n 
		query = oracle!(g, flow_val)
		if query == false
			return false
		end
		for e in g.edges
			e.wt *= (1 + (g.eps * abs(e.flow)) / (rho * e.cap))
			e.acc_flow += e.flow 
		end
	end
	for e in g.edges
		e.flow = (e.acc_flow * (1-g.eps)^2)/(n * (1+g.eps))
	end
	return true
end

# binary searches (with finite precision) on the max flow value which can be
# successfully approximated with a flow
# when it returns, the flow in the edges is the computed approx. max flow
# return type is void
function binary_search_max_flow!(g::Graph)
	low = g.b
	high = g.m*g.b
	# some part of g.eps is for error introduced in binary search over floats
	eps_bin = g.eps/2.0			# TODO: fine-tune
	# TODO: should g.eps be adjusted for constant in mult wt algo
	g.eps = g.eps - eps_bin		# TODO: actually g.eps should be some more for safety
	rho = 3 * sqrt(m/g.eps)
	while low + eps_bin*high < high
		mid = (low + high)/2.0
		query = mutl_wt!(g, mid, oracle!, rho)
		if query == false
			high = mid
		else
			low = mid
		end
	end
	mult_wt!(g, low, oracle!, rho)
end
