set.seed(105);
n = 200; p = 1000;

X = matrix(rnorm(n*p), nrow = n, ncol = p);
z1 = rnorm(n); z2 = rnorm(n);

X[, 1:4] = X[, 1:4]+z1;
X[, 5] = X[, 5]+2*z1;
X[, 6] = X[, 6]+1.5*z1;
X[, 7:20] = X[, 7:20]+0.5*z1;
X[, 21:40] = X[, 21:40]+0.5*z2;
beta = c(4, 2, -4, -2, rep(0, 996));
y = rnorm(n, X%*%beta, sd=1.5);

