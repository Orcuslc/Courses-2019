# library("hdrm");
# downloadData("Scheetz2006");
# attachData("Scheetz2006")

# (a)
for(alpha in c(1.0, 0.75, 0.5, 0.25)) {
  cvfit = cv.glmnet(X, y, family = "gaussian", alpha = alpha);
  lambda = cvfit$lambda.min;
  fit = glmnet(X, y, family = "gaussian", alpha = alpha, lambda = lambda);
  R_square = fit$dev.ratio;
  df = fit$df;
  print(c(alpha, R_square, df));
}

# (b)
for(alpha in c(1.0, 0.75, 0.5, 0.25)) {
  Chr5 = rownames(fData[which(fData["Chr"] == 5), ])
  index = match(Chr5, colnames(X))
  filtered_X = X[, index]
  
  cvfit = cv.glmnet(X, y, family = "gaussian", alpha = alpha);
  lambda = cvfit$lambda.min;
  fit = glmnet(X, y, family = "gaussian", alpha = alpha, lambda = lambda);
  R_square = fit$dev.ratio;
  df = fit$df;
  print(c(alpha, R_square, df));
}