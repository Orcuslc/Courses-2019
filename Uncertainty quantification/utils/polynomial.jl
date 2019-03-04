module Polynomial
export Legendre, Legendre_Normalized, Legendre_Normalized_Projection_Error, Legendre_Normalized_Projection_Coefficients, Legendre_Normalized_Projection, Monomial


function three_term_recurrence(x::Array, N::Integer, p0::Function, p1::Function, rec::Function)
	#=
		three term recurrence for orthogonal polynomials

		input:
		x: points, n*1 Array{Float64, 1}
		N: highest degree
		p0: first polynomial
		p1: second polynomial
		rec: recurrence relation rec(N, x, p0, p1)
	
		output:
		polys: n*(N+1) Array{Float64, 1}, containing polynomial values
	=#
	polys = zeros(size(x, 1), N+1);
	polys[:, 1] = p0(x);

	if N == 0
		return polys[:, 1];
	end

	polys[:, 2] = p1(x);
	for i in 2:N
	    polys[:, i+1] = rec(i-1, x, polys[:, i-1], polys[:, i]);
	end
	return polys;
end

function Monomial(x::Array, N::Integer)
	polys = zeros(size(x, 1), N+1);
	polys[:, 1] .= 1;
	for i in 1:N
		polys[:, i+1] = polys[:, i]	.* x;
	end
	return polys;
end

function Legendre(x::Array, N::Integer)
	p0 = x -> ones(size(x));
	p1 = x -> x
	rec = (n, x, p0, p1) -> (x.*p1*(2*n+1)-p0*n)/(n+1)
	return three_term_recurrence(x, N, p0, p1, rec)
end

function Legendre_Normalized(x::Array, N::Integer)
	p0 = x -> ones(size(x))*sqrt(1/2);
	p1 = x -> x*sqrt(3/2);
	rec = (n, x, p0, p1) -> (x.*p1*sqrt((2*n+1)*(2*n+3))-p0*n*sqrt((2*n+3)/(2*n-1)))/(n+1)
	return three_term_recurrence(x, N, p0, p1, rec)
end

function simpson(f::Function, interval::Array, n::Integer)
	#= 
		Composite Simpson quadrature
		f: target function
		interval: integral interval
		n: size of subintervals, should be even
	=#
	@assert(n%2 == 0);
	x = range(interval[1], stop = interval[2], length = n+1) |> collect;
	h = x[2]-x[1];
	return h/3*(f(x[1])+f(x[end])+4*sum(f.(x[2:2:end]))+2*sum(f.(x[3:2:end-1])));
end

function simpson(x::Array, h::Float64)
	#=
		Composite Simpson quadrature
		x: array of function values, the length of which should be odd
		h: step size
	=#
	@assert(size(x, 1) % 2 == 1);
	return h/3*(x[1]+x[end]+4*sum(x[2:2:end])+2*sum(x[3:2:end-1]));
end

function Legendre_Normalized_Projection_Coefficients(f::Function, N::Integer)
	#=
		Projection coefficients by Normalized Legendre polynomials. Computed through simpson quadrature.

		input:
		f: target function
		N: highest degree

		output:
		coeff: (N+1)*1 Array{Float64, 1}, coefficients of each polynomial
	=#
	n = 200000; # Large enough
	h = 2/n;
	x = range(-1., stop = 1., length = n+1) |> collect;
	legendre = Legendre_Normalized(x, N);
	func_values = legendre.*f.(x);
	sp = x -> simpson(x, h);
	coeff = mapslices(sp, func_values, dims = 1);
	return coeff[:];
end

function Legendre_Normalized_Projection(f::Function, x::Array, N::Integer)
	#=
		Approximation by Normalized Legendre polynomials.

		input:
		f: target function
		x: evaluation points
		N: highest degree

		output:
		y: size(x, 1)*1 Array{Float64, 1}, the approximation values;
	=#
	coeff = Legendre_Normalized_Projection_Coefficients(f, N);
	legendre = Legendre_Normalized(x, N);
	return sum(reshape(coeff, 1, :).*legendre, dims = 2);
end

function Legendre_Normalized_Projection_Error(f::Function, N::Integer)
	#=
		Approximation error of Normalized Legendre projection. Computed through composite simpson quadrature.
	=#
	n = 200000; # Large enough
	h = 2/n;
	x = range(-1., stop = 1., length = n+1) |> collect;
	func_values = (f.(x)-Legendre_Normalized_Projection(f, x, N)).^2;
	projection_error = simpson(func_values, h);
	return projection_error;
end

end # Polynomial