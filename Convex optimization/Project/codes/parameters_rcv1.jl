include("functions.jl")

# data
dataname = "rcv1.binary";
dataset = "../datasets/rcv1.binary/rcv1_sep.mat";

# number of hidden nodes
H = 128;

# number of iterations
K = 200;
checkpoint = 2;

# number of iterations for SSG
K_SSG = 1200;
checkpoint_SSG = 12;

# for SSG
batchsize = 2699;

# choise of functions
sigma = sigmoid;
dsigma = dsigmoid;
l = softmax;
dl = dsoftmax;

# for PG and APG, step size
x_eta = 4.0;
y_eta = 4.0;
z_eta = 4.0;
w_eta = 4.0;

# for APG
x_gamma0 = 1/x_eta;
y_gamma0 = 1/y_eta;
z_gamma0 = 1/z_eta;
w_gamma0 = 1/w_eta;

# for PG with line search
x_gamma_dec = 0.9;
x_gamma_inc = 0.5;
y_gamma_dec = 0.9;
y_gamma_inc = 0.5;
z_gamma_dec = 0.9;
z_gamma_inc = 0.5;
w_gamma_dec = 0.9;
w_gamma_inc = 0.5;

# for SSG, step size
C = 1.0;
