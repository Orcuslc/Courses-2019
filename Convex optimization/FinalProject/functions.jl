function relu(x)
	return max.(x, 0);
end

function sigmoid(x) 
	return exp.(x)./(1 .+ exp.(x));
end

function hinge(x)
	return max.(1 .- x, 0);
end 

function softmax(x)
	return log.(1 .+ exp.(-x));
end

# subgradient
function dsigmoid(x)
	return exp.(x) ./ (1 .+ exp.(x)).^2
end

function dhinge(x)
	if x > 0
		return 0 .* x;	# 0
	elseif x < 0
		return 0 .* x .- 1; # -1
	else
		return rand(Float64, size(x)) .* -1;
	end
end

function dsoftmax(x)
	return -exp.(-x) ./(1 .+ exp.(-x));
end