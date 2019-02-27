module Approximation
export legendre_discrete_projection

include("./polynomial.jl")
include("./quadrature.jl")
using .Polynomial:Legendre
using .Quadrature:legendre_gauss_quadrature

function legendre_discrete_projection(f::Function, x::Array, a::Float64, b::Float64, N::Integer)
	#=
		Discrete legendre projection I_N(f) = \sum_{i=1}^N \tilde{f_i} \phi_i(x), where \tilde{f_i} = 1/gamma_i*integral(f*phi_i*w).
		Input:
			f: target function
			x: computing nodes
			[a, b]: projection interval
			N: projection degree
		Output:
			pf: the projected function values
	=#
	gamma = map(n -> 2/(2*n+1), 0.:N) |> collect;
	phi = Legendre(x, N);
	print(phi)
	f_i = legendre_gauss_quadrature(x -> f(x).*Legendre(x, N), a, b, N+1);
	return phi*(1 ./gamma .* f_i);
end

function legendre_discrete_projection(f::Function, a::Float64, b::Float64, N::Integer)
	#=
		Discrete legendre projection I_N(f) = \sum_{i=1}^N \tilde{f_i} \phi_i(x), where \tilde{f_i} = 1/gamma_i*integral(f*phi_i*w).
		Input:
			f: target function
			[a, b]: projection interval
			N: projection degree
		Output:
			pf: the projected function
	=#
	gamma = map(n -> 2/(2*n+1), 0.:N) |> collect;
	phi = x -> Legendre(x, N);
	f_i = legendre_gauss_quadrature(x -> f(x).*Legendre(x, N), a, b, N+1);
	return x -> phi(x)*(1 ./gamma .* f_i);
end

end # module