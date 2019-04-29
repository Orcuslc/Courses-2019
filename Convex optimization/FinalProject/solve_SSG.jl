# include source code
include("neural_network.jl")
include("functions.jl")
include("parameters.jl")

using .parameters: dataname, dataset, H, K_SSG, checkpoint_SSG, batchsize, sigma, dsigma, l, dl;

# load data
using MAT
data = matread(dataset);
A = data["A_train"];
b = data["b_train"];
A_test = data["A_test"];
b_test = data["b_test"];

# proximal mapping; since the feasible set is R^n, it is just thd identity mapping
proj_x = x -> x;
proj_y = x -> x;
proj_z = x -> x;
proj_w = x -> x;

# sample xi
sample_xi = () -> rand(1:size(A, 1), (batchsize, 1))[:, 1];

# stepsize
C = 1.0;

# random initialize of parameters)
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

# first run
f, _, _ = forward(A, b, x, y, z, w, sigma, l);
output = predict(A_test, x, y, z, w, sigma);
acc = accuracy(output, b_test);

fpath[1, 1] = f;
accpath[1, 1] = acc;
println(f);
println(acc);

for i = 1:K_SSG

	# SSG
	eta = C/sqrt(i);
	index = sample_xi();

	global x, y, z, w;
	global dx, dy, dz, dw;
	global xsum, ysum, zsum, wsum;
	global xpath, ypath, zpath, wpath, dxpath, dypath, dzpath, dwpath, fpath, accpath;	
	
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

	# keep track of results
	if i % checkpoint_SSG == 0
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
h5write(string("./results/", dataname, "/SSG.h5"), "variables/x", xpath);
h5write(string("./results/", dataname, "/SSG.h5"), "variables/y", ypath);
h5write(string("./results/", dataname, "/SSG.h5"), "variables/z", zpath);
h5write(string("./results/", dataname, "/SSG.h5"), "variables/w", wpath);
h5write(string("./results/", dataname, "/SSG.h5"), "results/objective", fpath);
h5write(string("./results/", dataname, "/SSG.h5"), "results/dx", dxpath);
h5write(string("./results/", dataname, "/SSG.h5"), "results/dy", dypath);
h5write(string("./results/", dataname, "/SSG.h5"), "results/dz", dzpath);
h5write(string("./results/", dataname, "/SSG.h5"), "results/dw", dwpath);
h5write(string("./results/", dataname, "/SSG.h5"), "results/accuracy", accpath);