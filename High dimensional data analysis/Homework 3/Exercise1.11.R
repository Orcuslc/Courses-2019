load('~/GitHub/Spring-19-courses/High dimensional data analysis/Homework 3/whoari.RData');

variables = colnames(X);
form = as.formula(paste("y ~ ", paste(variables, collapse = "+")));

# normalization
normalize <- function(x) return((x-mean(x))/sd(x));
X = apply(X, 2, normalize);

# linear regression
lr_coeff = lm(form, as.data.frame(X));
sort(abs(lr_coeff$coefficients), decreasing = TRUE);

# ridge regression
ridge_regression <- function(X, y, lambda) {
	s = svd(X);
	U = s$u; V = s$v; D = s$d;
	rr_coeff = V %*% diag(D/(D^2+nrow(X)*lambda)) %*% t(U) %*% y;
	return(rr_coeff);
}

# generalized cross validation
GCV <- function(X, y, lambda) {
	s = svd(X);
	U = s$u; V = s$v; D = s$d;
	rr_coeff = V %*% diag(D/(D^2+nrow(X)*lambda)) %*% t(U) %*% y;
	beta0 = mean(y);
	RSS = sum((y - beta0 - X %*% rr_coeff)^2);
	df = sum(D^2/(D^2+nrow(X)*lambda));
	return(RSS/(1-1/nrow(X)*df)^2);
}

lambdas = seq(-5, 5, length=200);
lambdas = 10^lambdas;
GCVs = c();
for(i in 1:length(lambdas)) {
	GCV_ = GCV(X, y, lambdas[i]);
	GCVs = c(GCVs, GCV_);
}
index = sort(GCVs, index.return = TRUE)$ix[1];
lambda = lambdas[index];

# lambda = 0.1482021;
rr_coeff = ridge_regression(X, y, lambda);
names(rr_coeff) = variables;
variable_index = sort(abs(rr_coeff), decreasing = TRUE);
