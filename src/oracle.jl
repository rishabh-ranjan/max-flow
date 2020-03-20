include("struct_defs.jl")
using SparseArrays
using Laplacians

# (eps, 3*sqrt(m/eps)) oracle
# modifies the flows of the edges in-place
# returns false for "fail", else returns true
# Note: flows get corrupted (these are no more needed anyway)
# 		even when false is returned
function oracle!(g::Graph)

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
		e.flow = e.cond * (phi[e.to] - phi[e.from])
	end

	# calculate energy of the flow
	eng = 0.0
	for e in g.edges
		eng += e.flow^2 / e.cond
	end

	# TODO: move it inside the loop for micro-optimization
	return eng <= (1+g.eps)*ws
end
