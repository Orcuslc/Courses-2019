include("functions.jl")

# data
dataname = "bernoulli";
dataset = "bernoulli.mat";

# number of hidden nodes
H = 128;

# number of iterations
K = 200;
checkpoint = 2;

# number of iterations for SSG
K_SSG = 2000;
checkpoint_SSG = 20;

# for SSG
batchsize = 100;

# choise of functions
sigma = sigmoid;
dsigma = dsigmoid;
l = softmax;
dl = dsoftmax;

# for PG and APG, step size
x_eta = 0.4;
y_eta = 0.4;
z_eta = 0.4;
w_eta = 4.0;

# for APG
x_gamma0 = 1/(2*x_eta);
y_gamma0 = 1/(2*y_eta);
z_gamma0 = 1/(2*z_eta);
w_gamma0 = 1/(2*w_eta);

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
