# include source code
include("neural_network.jl")
include("functions.jl")

# load data
using MAT
data = matread("../datasets/rcv1.binary/rcv1.mat");
A = data["A"];
b = data["b"];

# proximal mapping; since the feasible set is R^n, it is just thd identity mapping
prox_x = (x, eta) -> x;
prox_y = (x, eta) -> x;
prox_z = (x, eta) -> x;
prox_w = (x, eta) -> x;

# parameters
sigma = sigmoid;
l = softmax;
dsigma = dsigmoid;
dl = dsoftmax;

# number of iterations
K = 100;

# random initialize of parameters
using Random
Random.seed!(12345);
H = 128;
x, y, z, w = initialize([size(A, 2), H, size(b, 2)]);

# solution path
xpath = Any[]; ypath = Any[]; zpath = Any[]; wpath = Any[]; fpath = Any[];
push!(xpath, x); push!(ypath, y); push!(zpath, z);push!(wpath, w);

# sum of solutions
xsum = x; ysum = y; zsum = z; wsum = w;

# derivative path
dxpath = Any[]; dypath = Any[]; dzpath = Any[]; dwpath = Any[];

# step size
x_eta = 4.0;
x_gamma0 = 1/x_eta;
x_z = x;
y_eta = 4.0;
y_gamma0 = 1/y_eta;
y_z = y;
z_eta = 4.0;
z_gamma0 = 1/z_eta;
z_z = z;
w_eta = 4.0;
w_gamma0 = 1/w_eta;
w_z = w;

# first run
f, _, _ = forward(A, b, x, y, z, w, sigma, l);
push!(fpath, f);
println(f)

for i = 1:K

	# using global variables in the local scope
	global x, y, z, w;
	global dx, dy, dz, dw;
	global xsum, ysum, zsum, wsum;	
	global x_gamma0, y_gamma0, z_gamma0, w_gamma0, x_eta, y_eta, z_eta, w_eta, x_z, y_z, z_z, w_z

	# compute the next step
	x_alpha = (-x_eta*x_gamma0 + sqrt((x_eta*x_gamma0)^2+4*x_eta*x_gamma0))/2;
	x_gamma1 = x_alpha^2/x_eta;
	x_y = 1/(x_alpha*x_gamma0+x_gamma1) .* (x_alpha*x_gamma0 .* x_z .+ x_gamma1 .* x);

	y_alpha = (-y_eta*y_gamma0 + sqrt((y_eta*y_gamma0)^2+4*y_eta*y_gamma0))/2;
	y_gamma1 = y_alpha^2/y_eta;
	y_y = 1/(y_alpha*y_gamma0+y_gamma1) .* (y_alpha*y_gamma0 .* y_z .+ y_gamma1 .* y);

	z_alpha = (-z_eta*z_gamma0 + sqrt((z_eta*z_gamma0)^2+4*z_eta*z_gamma0))/2;
	z_gamma1 = z_alpha^2/z_eta;
	z_y = 1/(z_alpha*z_gamma0+z_gamma1) .* (z_alpha*z_gamma0 .* z_z .+ z_gamma1 .* z);

	w_alpha = (-w_eta*w_gamma0 + sqrt((w_eta*w_gamma0)^2+4*w_eta*w_gamma0))/2;
	w_gamma1 = w_alpha^2/w_eta;
	w_y = 1/(w_alpha*w_gamma0+w_gamma1) .* (w_alpha*w_gamma0 .* w_z .+ w_gamma1 .* w);

	# forward pass at new values
	f, A1, A2 = forward(A, b, x_y, y_y, z_y, w_y, sigma, l);

	# backpropagation
	dx, dy, dz, dw = backward(f, A1, A2, A, b, x_y, y_y, z_y, w_y, sigma, l, dsigma, dl);

	# proximal mapping and update
	x1 = prox_x(x_y .- x_eta*dx, x_eta);
	x_z = x .+ 1/x_alpha*(x1 .- x);

	y1 = prox_y(y_y .- y_eta*dy, y_eta);
	y_z = y .+ 1/y_alpha*(y1 .- y);

	z1 = prox_z(z_y .- z_eta*dz, z_eta);
	z_z = z .+ 1/z_alpha*(z1 .- z);

	w1 = prox_w(w_y .- w_eta*dw, w_eta);
	w_z = w .+ 1/w_alpha*(w1 .- w);

	# update parameters
	x_gamma0 = x_gamma1;
	y_gamma0 = y_gamma1;
	z_gamma0 = z_gamma1;
	w_gamma0 = w_gamma1;

	x = x1;
	y = y1;
	z = z1;
	w = w1;

	# compute the average
	xsum += x1;
	ysum += y1;
	zsum += z1;
	wsum += w1;

	# save new result
	push!(xpath, xsum./(i+1));
	push!(ypath, ysum./(i+1));
	push!(zpath, zsum./(i+1));
	push!(wpath, wsum./(i+1));

	# compute the target
	f, _, _ = forward(A, b, xpath[end], ypath[end], zpath[end], wpath[end], sigma, l);
	push!(fpath, f);
	println(i);
	println(f);
en# include source code
include("neural_network.jl")
include("functions.jl")

# load data
using MAT
data = matread("../datasets/rcv1.binary/rcv1.mat");
A = data["A"];
b = data["b"];

# proximal mapping; since the feasible set is R^n, it is just thd identity mapping
prox_x = (x, eta) -> x;
prox_y = (x, eta) -> x;
prox_z = (x, eta) -> x;
prox_w = (x, eta) -> x;

# parameters
sigma = sigmoid;
l = softmax;
dsigma = dsigmoid;
dl = dsoftmax;

# number of iterations
K = 100;

# random initialize of parameters
using Random
Random.seed!(12345);
H = 128;
x, y, z, w = initialize([size(A, 2), H, size(b, 2)]);

# solution path
xpath = Any[]; ypath = Any[]; zpath = Any[]; wpath = Any[]; fpath = Any[];
push!(xpath, x); push!(ypath, y); push!(zpath, z);push!(wpath, w);

# sum of solutions
xsum = x; ysum = y; zsum = z; wsum = w;

# derivative path
dxpath = Any[]; dypath = Any[]; dzpath = Any[]; dwpath = Any[];

# step size
x_eta = 4.0;
x_gamma0 = 1/x_eta;
x_z = x;
y_eta = 4.0;
y_gamma0 = 1/y_eta;
y_z = y;
z_eta = 4.0;
z_gamma0 = 1/z_eta;
z_z = z;
w_eta = 4.0;
w_gamma0 = 1/w_eta;
w_z = w;

# first run
f, _, _ = forward(A, b, x, y, z, w, sigma, l);
push!(fpath, f);
println(f)

for i = 1:K

	# using global variables in the local scope
	global x, y, z, w;
	global dx, dy, dz, dw;
	global xsum, ysum, zsum, wsum;	
	global x_gamma0, y_gamma0, z_gamma0, w_gamma0, x_eta, y_eta, z_eta, w_eta, x_z, y_z, z_z, w_z

	# compute the next step
	x_alpha = (-x_eta*x_gamma0 + sqrt((x_eta*x_gamma0)^2+4*x_eta*x_gamma0))/2;
	x_gamma1 = x_alpha^2/x_eta;
	x_y = 1/(x_alpha*x_gamma0+x_gamma1) .* (x_alpha*x_gamma0 .* x_z .+ x_gamma1 .* x);

	y_alpha = (-y_eta*y_gamma0 + sqrt((y_eta*y_gamma0)^2+4*y_eta*y_gamma0))/2;
	y_gamma1 = y_alpha^2/y_eta;
	y_y = 1/(y_alpha*y_gamma0+y_gamma1) .* (y_alpha*y_gamma0 .* y_z .+ y_gamma1 .* y);

	z_alpha = (-z_eta*z_gamma0 + sqrt((z_eta*z_gamma0)^2+4*z_eta*z_gamma0))/2;
	z_gamma1 = z_alpha^2/z_eta;
	z_y = 1/(z_alpha*z_gamma0+z_gamma1) .* (z_alpha*z_gamma0 .* z_z .+ z_gamma1 .* z);

	w_alpha = (-w_eta*w_gamma0 + sqrt((w_eta*w_gamma0)^2+4*w_eta*w_gamma0))/2;
	w_gamma1 = w_alpha^2/w_eta;
	w_y = 1/(w_alpha*w_gamma0+w_gamma1) .* (w_alpha*w_gamma0 .* w_z .+ w_gamma1 .* w);

	# forward pass at new values
	f, A1, A2 = forward(A, b, x_y, y_y, z_y, w_y, sigma, l);

	# backpropagation
	dx, dy, dz, dw = backward(f, A1, A2, A, b, x_y, y_y, z_y, w_y, sigma, l, dsigma, dl);

	# proximal mapping and update
	x1 = prox_x(x_y .- x_eta*dx, x_eta);
	x_z = x .+ 1/x_alpha*(x1 .- x);

	y1 = prox_y(y_y .- y_eta*dy, y_eta);
	y_z = y .+ 1/y_alpha*(y1 .- y);

	z1 = prox_z(z_y .- z_eta*dz, z_eta);
	z_z = z .+ 1/z_alpha*(z1 .- z);

	w1 = prox_w(w_y .- w_eta*dw, w_eta);
	w_z = w .+ 1/w_alpha*(w1 .- w);

	# update parameters
	x_gamma0 = x_gamma1;
	y_gamma0 = y_gamma1;
	z_gamma0 = z_gamma1;
	w_gamma0 = w_gamma1;

	x = x1;
	y = y1;
	z = z1;
	w = w1;

	# compute the average
	xsum += x1;
	ysum += y1;
	zsum += z1;
	wsum += w1;

	# save new result
	push!(xpath, xsum./(i+1));
	push!(ypath, ysum./(i+1));
	push!(zpath, zsum./(i+1));
	push!(wpath, wsum./(i+1));

	# compute the target
	f, _, _ = forward(A, b, xpath[end], ypath[end], zpath[end], wpath[end], sigma, l);
	push!(fpath, f);
	println(i);
	println(f);
end