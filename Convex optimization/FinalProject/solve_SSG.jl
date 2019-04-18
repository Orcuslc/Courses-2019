# include source code
include("parameters.jl")
include("neural_network.jl")

# import packages
using .neural_network
using Statistics, LinearAlgebra

# load data
using MAT
data = matread("../datasets/news20/news20.mat")
A = data["X"]
b = data["y"]

# convert b to a sparse diagonal matrix
sparse_b = Diagonal(b[:, 1]);

# random initialize of parameters
x, y, z, w = initialize([size(A, 2), H, size(b, 2)]);

# proximal mapping; since the feasible set is R^n, it is just thd identity mapping
prox = x -> x;

# sample xi
sample_xi = () -> rand(1:size(x, 1), (batchsize, 1));

# functions
X1 = (x, y) -> A*x .+ transpose(y);
R = (x, y, z, w) -> b .* w + b .* sum(sigma(X1(x, y)) .* z, dims = 2);
f = (x, y, z, w) -> mean(l(R(x, y, z, w)), 1);