# include source code
include("neural_network.jl")
include("functions.jl")
include("parameters.jl")

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

x_z = x;
y_z = y;
z_z = z;
w_z = w;

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
	global x, y, z, w, f;
	global dx, dy, dz, dw;
	global xsum, ysum, zsum, wsum;	
	global x_gamma0, y_gamma0, z_gamma0, w_gamma0, x_eta, y_eta, z_eta, w_eta, x_z, y_z, z_z, w_z;
	global xpath, ypath, zpath, wpath, dxpath, dypath, dzpath, dwpath, fpath, accpath;

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
h5write(string("./results/", dataname, "/APG.h5"), "variables/x", xpath);
h5write(string("./results/", dataname, "/APG.h5"), "variables/y", ypath);
h5write(string("./results/", dataname, "/APG.h5"), "variables/z", zpath);
h5write(string("./results/", dataname, "/APG.h5"), "variables/w", wpath);
h5write(string("./results/", dataname, "/APG.h5"), "results/objective", fpath);
h5write(string("./results/", dataname, "/APG.h5"), "results/dx", dxpath);
h5write(string("./results/", dataname, "/APG.h5"), "results/dy", dypath);
h5write(string("./results/", dataname, "/APG.h5"), "results/dz", dzpath);
h5write(string("./results/", dataname, "/APG.h5"), "results/dw", dwpath);
h5write(string("./results/", dataname, "/APG.h5"), "results/accuracy", accpath);

# # draw convergence plot and accuracy plot
# using PyPlot;

# fig, ax1 = subplots()
# color = "tab:red"
# ax1.set_xlabel("iterations")
# ax1.set_ylabel('Objective', color = color)
# ax1.plot(fpath, color = color)
# ax1.tick_params(axis = 'y', labelcolor = color)

# ax2 = ax1.twinx()
# color = "tab:blue"
# ax2.set_ylabel('Accuracy (%)', color = color)
# ax2.plot(accpath, color = color)
# ax2.tick_params(axis = 'y', labelcolor = color)

# fig.tight_layout()
