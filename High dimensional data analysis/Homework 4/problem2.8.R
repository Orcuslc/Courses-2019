source("lasso.R");
source("cross_validation.R")
source("ridge_regression.R");

n = 50;

test <- function(X, y, beta) {
	lambdas = 10^seq(-2, 2, length = 20);
	lambda = lasso.cv(X, y, lambdas);
	beta_lasso = lasso.coordinate_descent(X, y, lambda);
	#beta_lasso = glmnet(X, y, lambda = lambda);
	error_lasso = sum((beta - beta_lasso)^2);
	
	lambdas = 10^seq(-2, 2, length=20);
	lambda = ridge_regression.cv(X, y, lambdas);
	beta_ridge = ridge_regression(X, y, lambda);
	error_ridge = sum((beta - beta_ridge)^2);
	
	X = as.data.frame(X);
	variables = colnames(X);
	form = as.formula(paste("y ~ 0+", paste(variables, collapse = "+")));
	fit = step(lm(form, data = X));
	beta_forward = as.data.frame(matrix(rep(0, length(beta)), nrow = 1));
	beta_forward[names(fit$coefficients)] = fit$coefficients;
	error_forward = sum((beta - beta_forward)^2);
	
	return(list(lasso=error_lasso, ridge=error_ridge, forward=error_forward));
}

# # problem 1
# forward_error = 0;
# lasso_error = 0;
# ridge_error = 0;
# for(i in 1:100) {
# 	p = 25;
# 	X = matrix(rnorm(n*p), nrow = n);
# 	beta = rep(0, p);
# 	beta[1] = 1; beta[2] = -1;
# 	y = rnorm(n, mean = X%*%beta);
# 	res = test(X, y, beta);
# 	forward_error = forward_error + res$forward;
# 	lasso_error = lasso_error + res$lasso;
# 	ridge_error = ridge_error + res$ridge;
# }
# forward_error = forward_error/100;
# lasso_error = lasso_error/100;
# ridge_error = ridge_error/100;
# print(c(forward_error, lasso_error, ridge_error))
# 
# # problem 2
# forward_error = 0;
# lasso_error = 0;
# ridge_error = 0;
# for(i in 1:100) {
# 	p = 100;
# 	X = matrix(rnorm(n*p), nrow = n);
# 	beta = rep(0, p);
# 	beta[1] = 1; beta[2] = -1;
# 	y = rnorm(n, mean = X%*%beta);
# 	res = test(X, y, beta);
# 	forward_error = forward_error + res$forward;
# 	lasso_error = lasso_error + res$lasso;
# 	ridge_error = ridge_error + res$ridge;
# }
# forward_error = forward_error/100;
# lasso_error = lasso_error/100;
# ridge_error = ridge_error/100;
# print(c(forward_error, lasso_error, ridge_error))

# # problem 3
# forward_error = 0;
# lasso_error = 0;
# ridge_error = 0;
# for(i in 1:100) {
# 	p = 25;
# 	X = matrix(rnorm(n*p), nrow = n);
# 	beta = rep(0, p);
# 	beta[1:4] = 0.5; beta[5:8] = -0.5;
# 	y = rnorm(n, mean = X%*%beta);
# 	res = test(X, y, beta);
# 	forward_error = forward_error + res$forward;
# 	lasso_error = lasso_error + res$lasso;
# 	ridge_error = ridge_error + res$ridge;
# }
# forward_error = forward_error/100;
# lasso_error = lasso_error/100;
# ridge_error = ridge_error/100;
# print(c(forward_error, lasso_error, ridge_error))

# problem 4
forward_error = 0;
lasso_error = 0;
ridge_error = 0;
for(i in 1:100) {
	p = 100;
	X = matrix(rnorm(n*p), nrow = n);
	beta = rep(0, p);
	beta[1:16] = 0.25; beta[17:32] = -0.25;
	y = rnorm(n, mean = X%*%beta);
	res = test(X, y, beta);
	forward_error = forward_error + res$forward;
	lasso_error = lasso_error + res$lasso;
	ridge_error = ridge_error + res$ridge;
}
forward_error = forward_error/100;
lasso_error = lasso_error/100;
ridge_error = ridge_error/100;
print(c(forward_error, lasso_error, ridge_error))