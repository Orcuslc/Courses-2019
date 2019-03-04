module Quadrature
export legendre_gauss_quadrature, legendre_gauss_weights_nodes

function legendre_gauss_weights_nodes(N::Integer, a::Float64, b::Float64)
	#=
		Legendre-Gauss quadrature weights and nodes
		Input:
			N: the number of quadrature points
			a, b: the quadrature interval
		Output:
			x: N*1 Array, the quadrature nodes
			w: N*1 Array, the quadrature weights
	=#
	# initial guess
	@assert(N >= 1);
	if(N == 1)
		x = [(b+a)/2];
		w = [(b-a)/2];
		return x, w;
	end
	x0 = cos.((2*collect(0.:1.:N-1).+1)*pi/(2*N))+(0.27/N)*sin.(pi.*range(-1, stop = 1., length = N).*(N-1)./(N+1)) |> collect;

	# Legendre-Gauss Vandermonde Matrix
	L = zeros(N, N+1);

	# devirative of L-G Vandermonde matrix (in fact, the last column)
	dL = zeros(N, 1);

	# Compute the zeros by recursion relation and Newton method
	x = zeros(N, 1) .+ 2;
	while(maximum(abs.(x0 .- x)) > eps())
		L[:, 1] .= 1;
		L[:, 2] .= x0;
		for k = 2:N
			L[:, k+1] .= ((2*k-1).*x0.*L[:, k] - (k-1).*L[:, k-1])./k;
		end
		dL = (N+1).*(L[:, N] .- x0.*L[:, N+1])./(1 .- x0.^2);
		x = x0;
		x0 = x - L[:, N+1]./dL;	
	end

	# map to [a, b]
	x = (a.*(1. .- x0) .+ b.*(1. .+ x0))./2;

	# weights
	w = (b-a)./((1 .- x0.^2).*dL.^2).*((N+1)/N)^2;
	return x, w;
end

function legendre_gauss_quadrature(f::Function, a::Float64, b::Float64, N::Integer)
	#=
		Legendre-Gauss quadrature
		Input:
			f: the target function
			[a, b]: quadrature interval
			N: the order of quadrature
		Output:
			I: the quadrature value
	=#
	x, w = legendre_gauss_weights_nodes(N, a, b);
	return sum(f(x).*w, dims = 1)';
end

end # module