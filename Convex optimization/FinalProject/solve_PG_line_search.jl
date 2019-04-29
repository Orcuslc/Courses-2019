# include source code
include("neural_network.jl")
include("functions.jl")
include("parameters.jl")

using .parameters: dataname, dataset, H, K, checkpoint, sigma, dsigma, l, dl;

# load data
using MAT
data = matread(dataset);
A = data["A_train"];
b = data["b_train"];
A_test = data["A_test"];
b_test = data["b_test"];

# proximal mapping; since the feasible set is R^n, it is just thd identity mapping
prox_x = (x, eta) -> x;
prox_y = (x, eta) -> x;
prox_z = (x, eta) -> x;
prox_w = (x, eta) -> x;

# random initialize of parameters
using Random
Random.seed!(12345);
x, y, z, w = initialize([size(A, 2), H, size(b, 2)]);

# solution path
xpath = x; ypath = y; zpath = z; wpath = w;

# sum of solutions
xsum = x; ysum = y; zsum = z; wsum = w;

# derivative path
dxpath = zeros(size(x)); dypath = zeros(size(y)); dzpath = zeros(size(z)); dwpath = zeros(size(w));

# objective and prediction accuracy path
fpath = zeros(1, 1);
accpath = zeros(1, 1);

# step size
x_eta = 4.0;
x_gamma_dec = 0.9;
x_gamma_inc = 0.5;
y_eta = 4.0;
y_gamma_dec = 0.9;
y_gamma_inc = 0.5;
z_eta = 4.0;
z_gamma_dec = 0.9;
z_gamma_inc = 0.5;
w_eta = 4.0;
w_gamma_dec = 0.9;
w_gamma_inc = 0.5;

# first run
f, _, _ = forward(A, b, x, y, z, w, sigma, l);
output = predict(A_test, x, y, z, w, sigma);
acc = accuracy(output, b_test);

fpath[1, 1] = f;
accpath[1, 1] = acc;
println(f);
println(acc);

for i = 1:K

	# using global variables in the local scope
	global x, y, z, w;
	global dx, dy, dz, dw;
	global xsum, ysum, zsum, wsum;	
	global x_eta, y_eta, z_eta, w_eta, x_gamma_dec, x_gamma_inc, y_gamma_inc, y_gamma_dec, z_gamma_inc, z_gamma_dec, w_gamma_inc, w_gamma_dec;
	global xpath, ypath, zpath, wpath, dxpath, dypath, dzpath, dwpath, fpath, accpath;

	# increase step size
	x_eta /= x_gamma_inc;
	y_eta /= y_gamma_inc;
	z_eta /= z_gamma_inc;
	w_eta /= w_gamma_inc;

	# forward pass
	f, A1, A2 = forward(A, b, x, y, z, w, sigma, l);

	# backward pass
	dx, dy, dz, dw = backward(f, A1, A2, A, b, x, y, z, w, sigma, l, dsigma, dl);

	# update x
	while true
		# set as global variable
		global x1;

		# decrease step size
		x_eta *= x_gamma_dec;

		# next step
		x1 = prox_x(x .- x_eta*dx, x_eta);

		# stop criterion, for x, regard it as a vector instead of a matrix
		f1, _, _ = forward(A, b, x1, y, z, w, sigma, l);
		x_vec = reshape(x, length(x));
		x1_vec = reshape(x1, length(x1));
		dx_vec = reshape(dx, length(dx));

		if f1 - f <= dx_vec' * (x1_vec - x_vec) + 1/(2*x_eta)*sum((x1_vec - x_vec).^2)
			break;
		end
	end

	# update y
	while true
		# set as global variable
		global y1;

		# decrease step size
		y_eta *= y_gamma_dec;

		# next step
		y1 = prox_y(y .- y_eta*dy, y_eta);

		# stop criterion
		f1, _, _ = forward(A, b, x, y1, z, w, sigma, l);
	    if f1 - f <= (dy' * (y1 - y))[1] + 1/(2*y_eta)*sum((y1 - y).^2)
	    	break;
	    end
	end

	# update z
	while true
		# set as global variable
		global z1;

		# decrease step size
		z_eta *= z_gamma_dec;

		# next step
		z1 = prox_z(z .- z_eta*dz, z_eta);

		# stop criterion
		f1, _, _ = forward(A, b, x, y, z1, w, sigma, l);
		if f1 - f <= (dz' * (z1 - z))[1] + 1/(2*z_eta)*sum((z1 - z).^2)
			break;
		end
	end

	# update w
	while true
		# set as global variable
		global w1;

		# decrease step size
		w_eta *= w_gamma_dec;

		# next step
		w1 = prox_w(w .- w_eta*dw, w_eta);

		# stop criterion
		f1, _, _ = forward(A, b, x, y, z, w1, sigma, l);
		if f1 - f <= (dw' * (w1 - w))[1] + 1/(2*w_eta)*sum((w1 - w).^2)
			break
		end
	end

	x = x1;
	y = y1;
	z = z1;
	w = w1;

	# compute the average
	xsum += x1;
	ysum += y1;
	zsum += z1;
	wsum += w1;

	# keep track of even iterations
	if i % checkpoint == 0
		# save new result
		xpath = cat(xpath, xsum./(i+1), dims = 3);
		ypath = cat(ypath, ysum./(i+1), dims = 3);
		zpath = cat(zpath, zsum./(i+1), dims = 3);
		wpath = cat(wpath, wsum./(i+1), dims = 3);

		# compute the objective value
		f, A1, A2 = forward(A, b, xpath[:, :, end], ypath[:, :, end], zpath[:, :, end], wpath[:, :, end], sigma, l);
		fpath = vcat(fpath, f);
		
		# compute the subgradient
		dx, dy, dz, dw = backward(f, A1, A2, A, b, xpath[:, :, end], ypath[:, :, end], zpath[:, :, end], wpath[:, :, end], sigma, l, dsigma, dl);
		dxpath = cat(dxpath, dx, dims = 3);
		dypath = cat(dypath, dy, dims = 3);
		dzpath = cat(dzpath, dz, dims = 3);
		dwpath = cat(dwpath, dw, dims = 3);

		# compute out-of-sample accuracy
		output = predict(A_test, xpath[:, :, end], ypath[:, :, end], zpath[:, :, end], wpath[:, :, end], sigma);
		acc = accuracy(output, b_test);
		accpath = vcat(accpath, acc);

		println(i);
		println(f);
		println(acc);
	end
end

# save all variables to file
using HDF5
using HDF5
h5write(string("./results/", dataname, "/PG1.h5"), "variables/x", xpath);
h5write(string("./results/", dataname, "/PG1.h5"), "variables/y", ypath);
h5write(string("./results/", dataname, "/PG1.h5"), "variables/z", zpath);
h5write(string("./results/", dataname, "/PG1.h5"), "variables/w", wpath);
h5write(string("./results/", dataname, "/PG1.h5"), "results/objective", fpath);
h5write(string("./results/", dataname, "/PG1.h5"), "results/dx", dxpath);
h5write(string("./results/", dataname, "/PG1.h5"), "results/dy", dypath);
h5write(string("./results/", dataname, "/PG1.h5"), "results/dz", dzpath);
h5write(string("./results/", dataname, "/PG1.h5"), "results/dw", dwpath);
h5write(string("./results/", dataname, "/PG1.h5"), "results/accuracy", accpath);