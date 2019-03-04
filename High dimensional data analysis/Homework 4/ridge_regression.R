# ridge regression
ridge_regression <- function(X, y, lambda) {
	s = svd(X);
	U = s$u; V = s$v; D = s$d;
	rr_coeff = V %*% diag(D/(D^2+nrow(X)*lambda)) %*% t(U) %*% y;
	return(rr_coeff);
}

GCV <- function(X, y, lambda) {
	s = svd(X);
	U = s$u; V = s$v; D = s$d;
	rr_coeff = V %*% diag(D/(D^2+nrow(X)*lambda)) %*% t(U) %*% y;
	beta0 = mean(y);
	RSS = sum((y - beta0 - X %*% rr_coeff)^2);
	df = sum(D^2/(D^2+nrow(X)*lambda));
	return(RSS/(1-1/nrow(X)*df)^2);
}

ridge_regression.cv <- function(X, y, lambdas) {
	GCVs = c();
	for(i in 1:length(lambdas)) {
		GCV_ = GCV(X, y, lambdas[i]);
		GCVs = c(GCVs, GCV_);
	}
	index = sort(GCVs, index.return = TRUE)$ix[1];
	lambda = lambdas[index];
	return(lambda);
}