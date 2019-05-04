library("hdrm")
library("splines")
library("grpreg")
downloadData("Golub1999")
attachData("Golub1999")

# group lasso
group_X = c();
for(i in 1:ncol(X)) {
  group_X = cbind(group_X, ns(X[, i], df = 3));
}
fit = cv.grpreg(group_X, y, rep(1:ncol(X), each = 3), family = "binomial");
lambda = which.min(fit$cve);
beta = fit$fit$beta[, lambda];
nparam = which(beta != 0)

# ordinary lasso
fit = cv.glmnet(X, y, family = "binomial");
index = which.min(fit$cvm);
beta = fit$glmnet.fit$beta[, index];
nparam = which(beta != 0)