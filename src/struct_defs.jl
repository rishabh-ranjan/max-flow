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
end
