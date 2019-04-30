include("neural_network.jl");
include("functions.jl");

using Random;

# number of hidden nodes
H = 128;

# dimension of each data point
d = 4;

# number of data samples
n = 1000;

# randomly choose x, y, z, w
Random.seed!(1000000007);
x, y, z, w = initialize([d, H, 1]);

# randomly generate A with a Gaussian distribution
Random.seed!(1234);
A = randn(Float64, (n, d));

# compute the probability of each b
sigma = sigmoid;
output = predict(A, x, y, z, w, sigma);
expoutput = exp.(output);
prob_1 = expoutput ./ (1 .+ expoutput);

# sample each b
# first, sample a random uniform t in (0, 1);
Random.seed!(123);
t = rand(Float64, (n, 1));

# second, if t_i > prob_i, then b_i = -1; otherwise b_i = 1;
function convert(x)
	if x == true
		return -1.0;
	else
		return 1.0;
	end
end

b = map(convert, t .> prob_1);

# save data to file
using MAT; 
matwrite("bernoulli.mat", Dict(
		"x" => x,
		"y" => y,
		"z" => z,
		"w" => w,
		"A_train" => A,
		"b_train" => b,
		"A_test" => A,
		"b_test" => b
	));