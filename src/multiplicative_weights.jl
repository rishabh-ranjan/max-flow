mutable struct Edge
	from		# src vertex id {1..n}
	to			# dest vertex id {1..n}
	cap			# capacity
	wt			# weight (used in mult wt)
	cond		# conductance (used in oracle)
	flow		# current flow (filled by oracle)
	acc_flow	# accumulated flow (used in mult wts)
end

function multiplicative_weights(s::Int,t::Int,n::Int,m::Int,flow,edge_list,oracle,eps,rho)
	for e in edge_list
		e.wt=1.0
		e.acc_flow=0.0
	end
	N=2*rho*(max(log(m),1))/(eps*eps)
	for i in 1:N 
		query=oracle(eps,n,m,s,t,flow,edge_list)
		if (query==false)
			return false
		end
		for e in edge_list:
			e.wt*=(1+(eps/rho)*(e.flow/e.cap))
			e.acc_flow+=e.flow 
			e.flow=0
		end
	end
	final_flow=0.0
	for e in edge_list:
		e.flow=((e.acc_flow*(1-eps)*(1-eps))/(N*(1+eps)))
		final_flow+=e.flow
	end
	return final_flow
end



function binary_search_max_flow(s::Int,t::Int,n::Int,m::Int,bottle_neck,edge_list,eps)
	low=bottle_neck
	high=m*bottle_neck
	eps1=eps/2.0 #change this!!
	eps2=eps-eps1
	rho=3*sqrt((m/eps2))
	while(low+eps1<high)
		mid=(low+high)/(2.0)
		query=multiplicative_weights(s,t,n,m,mid,edge_list,oracle,eps2,rho)
		if (query==false)
			high=mid
		else
			low=mid
		end
	end
	multiplicative_weights(s,t,n,m,low,edge_list,oracle,eps2,rho)
	return low
end