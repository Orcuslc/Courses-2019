lasso.coordinate_descent <- function(X, y, lambda) {
	# normalize X
	#normalize <- function(x) return((x-mean(x))/sd(x));
	#X = apply(X, 2, normalize);
	
	# number of variables
	p = ncol(X);
	
	# initialize beta with 0
	beta = rep(0, p);
	
	# compute residue
	r = y - X %*% beta;
	
	# tolerance
	tol = 1e-8;
	
	# iteration
	beta1 = rep(1, p);
	while(1) {
		for(j in 1:p) {
			# OLS estimator
			z = sum(X[, j]*r)/nrow(X) + beta[j];
			
			# S(z|lambda)
			if(z > lambda) {
				beta1[j] = z - lambda;
			}
			else if(z < -lambda) {
				beta1[j] = z + lambda;
			}
			else {
				beta1[j] = 0;
			}
			
			# update residue
			r = r - (beta1[j] - beta[j]) * X[, j];
		}
		if(sqrt(sum((beta - beta1)^2)) < tol)
			break;
		
		# update beta
		beta = beta1;
	}
	return(beta1);
}