module parameters
export dataname, dataset, H, K, batchsize, sigma, dsigma, l, dl;

include("functions.jl")

# data
# dataname = "rcv1.binary";
# dataset = "../datasets/rcv1.binary/rcv1_sep.mat";

dataname = "covtype";
dataset = "../datasets/covtype/covtype_sep.mat";

# number of hidden nodes
H = 128;

# number of iterations
K = 200;
checkpoint = 2;

# number of iterations for SSG
K_SSG = 2000;
checkpoint_SSG = 200;

# for SSG
batchsize = 46481;

# choise of functions
sigma = sigmoid;
dsigma = dsigmoid;
l = softmax;
dl = dsoftmax;

end # module