include("../utils/approximation.jl");
include("../utils/quadrature.jl");
using .Approximation;
using .Quadrature:legendre_gauss_weights_nodes
using LinearAlgebra;
using Random;
Random.seed!(1);

x0 = range(-1., stop = 1., length = 1000) |> collect;

f1 = x -> abs.(x);
f2 = x -> abs.(cos.(pi.*x)).^3;
f3 = x -> cos.(pi.*x);
f4 = x -> 1. ./(1. .+ 25. .*x.^2);
funcs = [f1, f2, f3, f4];

# # Problem 1, discrete projection
# err1 = zeros(21, 4)
# for i = 1:21
# 	for k = 1:4
# 		err1[i, k] = norm(legendre_discrete_projection(funcs[k], -1., 1., i-1)(x0) - funcs[k](x0));
# 	end
# end

# # Problem 2, Legendre interpolation with legendre-gauss quadrature points
err2 = zeros(21, 4)
for i = 1:21
	x, w = legendre_gauss_weights_nodes(i, -1., 1.);
	print(size(x))
	for k = 1:4
		err2[i, k] = norm(legendre_interpolation(x, funcs[k](x))(x0) - funcs[k](x0));
	end
end

# # Problem 3, Legendre interpolation with equispaced points
# err3 = zeros(21, 4)
# for i = 1:21
# 	if(i == 1)
# 		x = [0.];
# 	else
# 		x = range(-1., stop = 1., length = i) |> collect;
# 	end
# 	for k = 1:4
# 		err3[i, k] = norm(lagrange_interpolation(x, funcs[k](x))(x0) - funcs[k](x0));
# 	end
# end

# # Problem 4, LS approximation with uniform random nodes
# err4 = zeros(21, 4)
# for i = 1:21
# 	x = rand(2*i, 1).*2 .- 1.;
# 	for k = 1:4
# 		err4[i, k] = norm(legendre_least_square(x, funcs[k](x), i-1)(x0) - funcs[k](x0));
# 	end
# end

# Problem 5, LS approximation with uniform random nodes
# err5 = zeros(21, 4)
# for i = 1:21
# 	x = rand(2*i, 1).*2 .- 1.;
# 	for k = 1:4
# 		err5[i, k] = norm(monomial_least_square(x, funcs[k](x), i-1)(x0) - funcs[k](x0));
# 	end
# end

using PyCall;
@pyimport matplotlib.pyplot as plt
labels = ["|x|", "|cos(pi*x)|^3", "cos(pi*x)", "1/(1+25x^2)"]

# # Problem 1
# N = 0:20;
# plt.figure();
# for k = 1:4
# 	plt.semilogy(N, err1[:, k], label = labels[k]);
# end
# plt.legend();
# plt.title("discrete legendre projection")
# plt.xlabel("N")
# plt.ylabel("error")
# plt.grid(true)
# plt.show()

# # Problem 2
N = 0:20;
plt.figure();
for k = 1:4
	plt.semilogy(N, err2[:, k], label = labels[k]);
end
plt.legend();
plt.title("legendre interpolation with quadrature points")
plt.xlabel("N")
plt.ylabel("error")
plt.grid(true)
plt.show()

# # Problem 3
# N = 0:20;
# plt.figure();
# for k = 1:4
# 	plt.semilogy(N, err3[:, k], label = labels[k]);
# end
# plt.legend();
# plt.title("legendre interpolation with uniform points")
# plt.xlabel("N")
# plt.ylabel("error")
# plt.grid(true)
# plt.show()

# # Problem 4
# N = 0:20;
# plt.figure();
# for k = 1:4
# 	plt.semilogy(N, err4[:, k], label = labels[k]);
# end
# plt.legend();
# plt.title("LS approximation with uniform random nodes")
# plt.xlabel("N")
# plt.ylabel("error")
# plt.grid(true)
# plt.show()

# # Problem 5
# N = 0:20;
# plt.figure();
# for k = 1:4
# 	plt.semilogy(N, err5[:, k], label = labels[k]);
# end
# plt.legend();
# plt.title("LS approximation with uniform random nodes and monomial basis")
# plt.xlabel("N")
# plt.ylabel("error")
# plt.grid(true)
# plt.show()
