module optimizer

function proximal_gradient(Dg, prox, x0, eta, N, epsilon)
	#= 
		Proximal Gradient without line search for f = g+h;
		Input:
			Dg: gradient of g
			prox: the proximal mapping w.r.t. eta*h
			x0: initial guess
			eta: step size
			epsilon: terminate tolerance
	=#
	x_path = x0;
	while(1)
		x1 = prox(x0 - eta*Dg(x0));
	end	
end

function stochastic_subgradient(G, prox, C, sample_xi, x0, K)
	#=
		Stochastic subgradient method 
		Input:
			G: stochastic subgradient function G(x, xi)
			prox: proximal mapping
			C: constant for control stepsize
			sample_xi: the function to sample xi
			x0: initial guess
			K: number of iterations
		Output:
			xpath: iteration path of average values of x
	=#
	xpath = x0;
	xsum = x0;
	for i = 1:K
		eta = C/sqrt(i);
		xi = sample_xi();
		x0 = prox(x0 - eta*G(x0, xi), eta);
		xsum = xsum + x0;
		xpath = hcat(xpath, xsum/(i+1));
	end
	return xpath;
end
	
end # module