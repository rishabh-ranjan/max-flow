using SparseArrays
using Laplacians

mutable struct Edge
	from		# src vertex id {1..n}
	to			# dest vertex id {1..n}
	cap			# capacity
	wt			# weight (used in mult wt)
	cond		# conductance (used in oracle)
	flow		# current flow (filled by oracle)
	acc_flow	# accumulated flow (used in mult wts)
end

# (eps, 3*sqrt(m/eps)) oracle
# modifies the flows of the edges in-place
# returns false for "fail", else returns true
# Note: flows get corrupted (these are no more needed anyway)
# 		even when false is returned
# eps -> global epsilon
# n -> number of vertices
# m -> number of edges
# s -> source vertex
# t -> sink vertex
# edges -> array of edges 
function oracle!(eps, n, m, s, t, edges::Array{Edge, 1})

	# sum of wts
	ws = 0.0
	for e in edges
		ws += e.wt
	end
	# calculation of conductances
	for e in edges
		e.cond = e.cap^2 / (e.wt + (eps*ws)/(3*m))
	end

	# construction of laplacian matrix
	lap = SparseArrays.spzeros(n, n)
	for e in edges
		lap[e.from, e.to] = -e.cond
		lap[e.to, e.from] = -e.cond
		lap[e.from, e.from] += e.cond
		lap[e.to, e.to] += e.cond
	end

	# construction of chi_s,t
	# TODO: this need not be recomputed, move it out
	chi = SparseArrays.spzeros(n)
	chi[s] = 1
	chi[t] = -1

	# TODO: pass the del = eps/3 parameter
	solver = Laplacians.KMPSDDMSolver(lap)
	# vertex potentials
	phi = solver(chi)

	for e in edges
		e.flow = e.cond * (phi[e.to] - phi[e.from])
	end

	# calculate energy of the flow
	eng = 0.0
	for e in edges
		eng += e.flow^2 / e.cond
	end

	# TODO: move it inside the loop for micro-optimization
	return eng <= (1+eps)*ws
end
