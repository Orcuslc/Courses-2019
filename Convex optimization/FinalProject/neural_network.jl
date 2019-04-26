include("functions.jl")
include("optimizers.jl")

function initialize(dimensions)
	#=
		Initialize parameters for the network
		Input:
			dimensions: Array, [input, hidden, output (=1 here)]
	=#
	x = rand(Float64, (dimensions[1], dimensions[2])) .* 2 .- 1;
	y = rand(Float64, dimensions[2]) .* 2 .- 1;
	z = rand(Float64, (dimensions[2], dimensions[3])) .*2 .- 1;
	w = rand(Float64, dimensions[3]) .*2 .- 1;
	return x, y, z, w;
end

function forward(A, b, x, y, z, w, sigma, l)
	# two-layer neural network
	# A: data matrix
	# b: label vector, represented as sparse diagonal matrix here
	# x: weight for the first layer
	# y: bias for the first layer
	# z: wieght for the second layer
	# w: bias for the second layer
	# sigma: activation function, relu or sigmoid
	# l: output function, hinge or softmax
	A1 = A*x .+ y';
	A2 = sigma.(A1)*z .+ w;
	f = sum(l.(b .* A2), dims = 1)/size(b, 1);
	return f[1], A1, A2;
end

function backward(f, A1, A2, A, b, x, y, z, w, sigma, l, dsigma, dl)
	#=
		Get the gradients for each variables by backpropagation
		Input:
			f: loss, defined in (2.2)
			A1, A2: output of each layer, defined in (2.2)
			A, b: training data and label
			x, y, z, w: network weights and biases
			sigma, l: activation function and loss function
			dsigma, dl: derivative of sigma and l
	=#

	# number of samples
	n = size(b, 1);

	# derivative for the second layer
	dA2 = dl.(b .* A2) .* b ./ n;
	dw = dA2'*ones(n, 1);
	dz = sigma.(A1)' * dA2;

	# derivative for the first layer
	dA1 = dA2*z' .* dsigma.(A1);
	dy = dA1'*ones(n, 1);
	dx = A'*dA1;
	return dx, dy, dz, dw;
end

# function compute_R(A, x, y)
# 	#=
# 		defined in (2.4)
# 	=#
# 	return A*x .+ y';
# end

# function compute_W(R, b, z, w, sigma)
# 	#= 
# 		defined in (2.4)

# 		Input:
# 			- R: defined in (2.4), computed by `compute_R`
# 			- b: sparse diagonal matrix
# 	=#
# 	return b*(w .+ sigma(R)*z);
# end

# function dw(W, b, dl)
# 	#=
# 		subgradient of w, defined in (2.6) 

# 		Input:
# 			- W: defined in (2.4), computed by `forward`
# 			- b: column vector, not sparse matrix
# 			- dl: subgradient of l
# 	=#
# 	return b'*dl(W) / size(b, 1);
# end

# function dz(R, W, b, sigma, dl)
	 
# 		subgradient of z, defined in (2.6)

# 		Input:
# 			- R, W: defined in (2.4), computed by `forward`
# 			- b: column vector, not sparse matrix
# 			- dl: subgradient of l
	
# 	return sigma(R)'*(dl(W) .* b) / size(b, 1);
# end

# function dy(R, W, b, z, dsigma, dl)
# 	#= 
# 		subgradient of z, defined in (2.6)

# 		Input:
# 			- R, W: defined in (2.4), computed by `forward`
# 			- b: column vector, not sparse matrix
# 			- dsigma, dl: subgradient of sigma, l
# 	=#
# 	return dsigma(R)'*(dl(W) .* b) .* z / size(b, 1);
# end

# function dx(R, W, A, b, z, dsigma, dl)
# 	#=
# 		subgradient of x, defined in (2.6)

# 		Input:
# 			- R, W: defined in (2.4), computed by `forward`
# 			- b: column vector, not sparse matrix
# 			- dsigma, dl: subgradient of sigma, l
# 	=#
# 	tmp = dl(W) .* b .* dsigma(R);
# 	return (Diagonal(tmp[:, 1])*A)' .* z' / size(b, 1);
# end
