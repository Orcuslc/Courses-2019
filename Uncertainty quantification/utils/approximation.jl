module Approximation
export legendre_discrete_projection, lagrange_interpolation

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

function newton_interpolation_coefficients(x::Array, fx::Array)
	#= 
		Interpolation coefficients for newton interpolation scheme
		Input:
			x: interpolation nodes
			fx: function values on corresponding node
		Output:
			c: coefficient for each interpolation basis function
	=#
	N = size(x, 1);
	for i = 1:(N-1)
		fx[i+1:end] = (fx[i+1:end] .- fx[i:end-1]) ./ (x[i+1:end] .- x[1:end-i]);
	end	
	return fx;
end

function horner_scheme(x::Array, c::Array, x0::Array)
	#=
		Horner scheme for evaluating a polynomial with given x and coefficient
		Input:
			x: evaluation nodes
			c: polynomial coefficients, [c0, c1, ..., cN]
			x0: given interpolation nodes
	=#
	y = c[end];
	for i = length(c)-1:-1:1
		y = y.*(x .- x0[i]) .+ c[i];
	end
	return y;
end

function lagrange_interpolation(x::Array, fx::Array)
	#=
		lagrange interpolation via Newton interpolation formula
		Input:
			x: interpolation nodes, n*1 Array{Float64, 1}
			fx: corresponding function values, n*1 Array{Float64, 1}
		Output:
			Pif: the function of degree-(n-1) lagrange interpolation
	=#
	c = newton_interpolation_coefficients(x, fx);
	return x_eval -> horner_scheme(x_eval, c, x);
end

function lagrange_interpolation(f::Function, x::Array)
#=
		lagrange interpolation via Newton interpolation formula
		Input:
			f: the target function
			x: interpolation nodes, n*1 Array{Float64, 1}
		Output:
			Pif: the function of degree-(n-1) lagrange interpolation
	=#
	c = newton_interpolation_coefficients(x, f(x));
	return x_eval -> horner_scheme(x_eval, c, x);
end

end # module