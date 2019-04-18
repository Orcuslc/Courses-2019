include("functions.jl")
using .functions

# number of hidden nodes
H = 128;

# number of iterations
K = 1000;

# for SSG
batchsize = 1; 

# choise of functions
l = hinge;
dl = dhinge;