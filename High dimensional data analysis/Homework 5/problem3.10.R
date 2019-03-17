library(splines);
attachData(Koussounadis2014);
sDay <- ns(sData$Day, df=2);
X0 <- model.matrix(~ sData$Treatment*sDay, sData)[,-1];
w <- rep(0:1, c(ncol(X0), ncol(X)));
XX <- cbind(X0, X)

# gamma = seq(2.1, 4.0, 0.1);
# err = c();
# for(gam in gamma) {
#   fit = cv.ncvreg(XX, y, penalty = "MCP", gamma = gam);
#   err = c(err, min(fit$cve))
# }
# gamma = gamma[which.min(err)];
# print(gamma)

# fit = cv.ncvreg(XX, y, penalty = "MCP", gamma = 2.3);
# plot(fit, type = 'rsq');
# 
# fit = cv.ncvreg(XX, y, penalty = "lasso");
# plot(fit, type = 'rsq');

fit = ncvreg(XX, y, penalty = "MCP", gamma = 2.3) # 2038
plot(fit)

fit = ncvreg(XX, y, penalty = "lasso") # 7912
plot(fit)

for(n in 1:100) {
  for(g in 10:nrow(fit$beta)) {
    if(fit$beta[g, n] != 0) {
      print(c(n, g));
      break;
    }
  }
}

lines(fit$lambda, fit$beta[2038,], xlim = c(max(fit$lambda), min(fit$lambda)));
lines(fit$lambda, fit$beta[7912,], xlim = c(max(fit$lambda), min(fit$lambda)))