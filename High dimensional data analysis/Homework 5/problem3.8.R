p <- function(theta, lambda, tao) {
  return(lambda^2/tao*(1-exp(-tao*abs(theta)/lambda)));
}

dp <- function(theta, lambda, tao) {
  return(exp(-tao/lambda*abs(theta)));
}

lambda = 2;
tao = 0.3;

theta = seq(-lambda/tao, lambda/tao, 0.01);
penalty = p(theta, lambda, tao);
plot(theta, penalty, xaxt="n", yaxt="n");

theta = seq(0, lambda/tao, 0.01);
Dp = dp(theta, lambda, tao);
plot(theta, Dp, xaxt="n", yaxt="n");
