n = 50;

test <- function(X, y, beta) {
  # lasso
  fit = cv.ncvreg(X, y);
  beta_lasso = fit$fit$beta[, which.min(fit$cve)]
  error_lasso = sum((beta - beta_lasso[2:length(beta_lasso)])^2);
  
  # ridge
  library("hdrm");
  fit = ridge(X, y);
  beta_ridge = fit$beta[, which.min(fit$GCV)];
  error_ridge = sum((beta - beta_ridge[2:length(beta_ridge)])^2);
  
  # forward
  X = as.data.frame(X);
  variables = colnames(X);
  form = as.formula(paste("y ~ 0+", paste(variables, collapse = "+")));
  fit = step(lm(form, data = X), k = log(nrow(X)));
  beta_forward = as.data.frame(matrix(rep(0, length(beta)), nrow = 1));
  beta_forward[names(fit$coefficients)] = fit$coefficients;
  error_forward = sum((beta - beta_forward)^2);
  
  # MCP
  fit = cv.ncvreg(X, y, penalty = "MCP");
  beta_MCP = fit$fit$beta[, which.min(fit$cve)];
  error_MCP = sum((beta - beta_MCP[2:length(beta_MCP)])^2);
  
  # SCAD
  fit = cv.ncvreg(X, y, penalty = "SCAD");
  beta_SCAD = fit$fit$beta[, which.min(fit$cve)];
  error_SCAD = sum((beta - beta_SCAD[2:length(beta_SCAD)])^2);
  
  return(list(lasso=error_lasso, ridge=error_ridge, forward=error_forward, MCP=error_MCP, SCAD=error_SCAD));
}

test_without_forward <- function(X, y, beta) {
  # lasso
  fit = cv.ncvreg(X, y);
  beta_lasso = fit$fit$beta[, which.min(fit$cve)]
  error_lasso = sum((beta - beta_lasso[2:length(beta_lasso)])^2);
  
  # ridge
  library("hdrm");
  fit = ridge(X, y);
  beta_ridge = fit$beta[, which.min(fit$GCV)];
  error_ridge = sum((beta - beta_ridge[2:length(beta_ridge)])^2);

  # MCP
  fit = cv.ncvreg(X, y, penalty = "MCP");
  beta_MCP = fit$fit$beta[, which.min(fit$cve)];
  error_MCP = sum((beta - beta_MCP[2:length(beta_MCP)])^2);
  
  # SCAD
  fit = cv.ncvreg(X, y, penalty = "SCAD");
  beta_SCAD = fit$fit$beta[, which.min(fit$cve)];
  error_SCAD = sum((beta - beta_SCAD[2:length(beta_SCAD)])^2);
  
  return(list(lasso=error_lasso, ridge=error_ridge, MCP=error_MCP, SCAD=error_SCAD));
}

# problem 1
forward_error = 0;
lasso_error = 0;
ridge_error = 0;
mcp_error = 0;
scad_error = 0;
for(i in 1:100) {
	p = 25;
	X = matrix(rnorm(n*p), nrow = n);
	beta = rep(0, p);
	beta[1] = 1; beta[2] = -1;
	y = rnorm(n, mean = X%*%beta);
	res = test(X, y, beta);
	forward_error = forward_error + res$forward;
	lasso_error = lasso_error + res$lasso;
	ridge_error = ridge_error + res$ridge;
	mcp_error = mcp_error + res$MCP;
	scad_error = scad_error + res$SCAD;
}
print(c(forward_error/100, lasso_error/100, ridge_error/100, mcp_error/100, scad_error/100))


# problem 2
forward_error = 0;
lasso_error = 0;
ridge_error = 0;
mcp_error = 0;
scad_error = 0;
for(i in 1:100) {
	p = 100;
	X = matrix(rnorm(n*p), nrow = n);
	beta = rep(0, p);
	beta[1] = 1; beta[2] = -1;
	y = rnorm(n, mean = X%*%beta);
	res = test_without_forward(X, y, beta);
	lasso_error = lasso_error + res$lasso;
	ridge_error = ridge_error + res$ridge;
	mcp_error = mcp_error + res$MCP;
	scad_error = scad_error + res$SCAD;
}
print(c(forward_error/100, lasso_error/100, ridge_error/100, mcp_error/100, scad_error/100))


# problem 3
forward_error = 0;
lasso_error = 0;
ridge_error = 0;
mcp_error = 0;
scad_error = 0;
for(i in 1:100) {
	p = 25;
	X = matrix(rnorm(n*p), nrow = n);
	beta = rep(0, p);
	beta[1:4] = 0.5; beta[5:8] = -0.5;
	y = rnorm(n, mean = X%*%beta);
	res = test(X, y, beta);
	forward_error = forward_error + res$forward;
	lasso_error = lasso_error + res$lasso;
	ridge_error = ridge_error + res$ridge;
	mcp_error = mcp_error + res$MCP;
	scad_error = scad_error + res$SCAD;
}
print(c(forward_error/100, lasso_error/100, ridge_error/100, mcp_error/100, scad_error/100))

# problem 4
forward_error = 0;
lasso_error = 0;
ridge_error = 0;
mcp_error = 0;
scad_error = 0;
for(i in 1:100) {
  p = 100;
  X = matrix(rnorm(n*p), nrow = n);
  beta = rep(0, p);
  beta[1:16] = 0.25; beta[17:32] = -0.25;
  y = rnorm(n, mean = X%*%beta);
  res = test_without_forward(X, y, beta);
  forward_error = forward_error + res$forward;
  lasso_error = lasso_error + res$lasso;
  ridge_error = ridge_error + res$ridge;
  mcp_error = mcp_error + res$MCP;
  scad_error = scad_error + res$SCAD;
}
print(c(forward_error/100, lasso_error/100, ridge_error/100, mcp_error/100, scad_error/100))