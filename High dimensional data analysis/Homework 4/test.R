## test coordinate

source("lasso.R");
pollution <- readRDS("~/GitHub/Spring-19-courses/High dimensional data analysis/Homework 4/pollution.rds");
X = pollution$X;
y = pollution$y;

# normalize X
normalize <- function(x) return((x-mean(x))/sd(x));
X = apply(X, 2, normalize);

# lambdas for test
lambdas = 10^seq(-2, 2, length = 10);

# test
for(i in 1:length(lambdas)) {
	beta_lasso = lasso.coordinate_descent(X, y, lambdas[i]);
	beta_glm = glmnet(X, y, lambda = lambdas[i])$beta;
	print(sum((beta_glm - beta_lasso)^2));
}




## test cross validation

source("cross_validation.R");
source("lasso.R");

load("~/GitHub/Spring-19-courses/High dimensional data analysis/Homework 4/whoari.RData");

# normalize X
normalize <- function(x) return((x-mean(x))/sd(x));
X = apply(X, 2, normalize);

lambdas = 10^seq(-2, 0, length = 20);
lambda = lasso.cv(X, y, lambdas);
print(lambda)
beta = lasso.coordinate_descent(X, y, lambda)

## check lasso is better ?

X_train = X[1:500, ]
y_train = y[1:500]
X_test = X[501:816, ]
y_test = y[501:816]
lambdas = 10^seq(-2, 0, length = 20);
lambda = lasso.cv(X_train, y_train, lambdas);
print(lambda)
beta_lasso = lasso.coordinate_descent(X_train, y_train, lambda);
s = svd(X_train);
beta_OLS = s$v %*% diag(1/s$d) %*% t(s$u) %*% y_train;
beta_null = rep(0, length(beta_OLS));

MSPE_lasso = mean((y_test - X_test%*%beta_lasso)^2);
MSPE_OLS = mean((y_test - X_test%*%beta_OLS)^2);
MSPE_null = mean((y_test - X_test%*%beta_null)^2);
print(c(MSPE_lasso, MSPE_OLS, MSPE_null));
