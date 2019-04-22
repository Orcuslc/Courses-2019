# include source code
include("neural_network.jl")
include("functions.jl")

# load data
using MAT
data = matread("../datasets/rcv1.binary/rcv1.mat");
A = data["A"];
b = data["b"];

# proximal mapping; since the feasible set is R^n, it is just thd identity mapping
proj_x = x -> x;
proj_y = x -> x;
proj_z = x -> x;
proj_w = x -> x;

# sample xi
batchsize = 100;
sample_xi = () -> rand(1:size(A, 1), (batchsize, 1))[:, 1];

# parameters
sigma = sigmoid;
l = softmax;
dsigma = dsigmoid;
dl = dsoftmax;

# number of iterations
K = 100;

# stepsize
C = 1.0;

# random initialize of parameters
H = 128;
x, y, z, w = initialize([size(A, 2), H, size(b, 2)]);

# solution path
xpath = Any[]; ypath = Any[]; zpath = Any[]; wpath = Any[]; fpath = Any[];
push!(xpath, x); push!(ypath, y); push!(zpath, z);push!(wpath, w);

# sum of solutions
xsum = x; ysum = y; zsum = z; wsum = w;

# derivative path
dxpath = Any[]; dypath = Any[]; dzpath = Any[]; dwpath = Any[];

for i = 1:K

	# SSG
	eta = C/sqrt(i);
	index = sample_xi();

	global x, y, z, w;
	global dx, dy, dz, dw;
	global xsum, ysum, zsum, wsum;
	
	# forward pass
	f, A1, A2 = forward(A[index, :], b[index, :], x, y, z, w, sigma, l);

	# backpropagation
	dx, dy, dz, dw = backward(f, A1, A2, A[index, :], b[index, :], x, y, z, w, sigma, l, dsigma, dl);

	# update
	x = proj_x(x - eta*dx);
	y = proj_y(y - eta*dy);
	z = proj_z(z - eta*dz);
	w = proj_w(w - eta*dw);

	# compute the average
	xsum += x;
	ysum += y;
	zsum += z;
	wsum += w;

	# save new result
	push!(xpath, x./(i+1));
	push!(ypath, y./(i+1));
	push!(zpath, z./(i+1));
	push!(wpath, w./(i+1));

	# compute the target
	f, _, _ = forward(A, b, xpath[end], ypath[end], zpath[end], wpath[end], sigma, l);
	push!(fpath, f);
	println(i);
	println(f);
end