ex1.2 <- function (N = 1000, n = 50, p = 5) {
	set.seed(1e9+7);
	# outputs
	res <- mse <- mspe <- coverage <- NULL;
	conf_int <- NULL;
	# progress bar
	pb <- txtProgressBar(0, N, style = 3);
	for(i in 1:N) {
		# for fit
		X <- matrix(runif(n * p), n, p);
		y <- rnorm(n);
		
		# for predict
		XX <- matrix(runif(n * p), n, p);
		yy <- rnorm(n);
		
		# marginal regression
		coeff <- NULL;
		for(j in 1:p) {
			vData = as.data.frame(matrix(X[, j]));
			fit = lm(y ~ V1+1, data = vData);
			coeff = c(coeff, summary(fit)$coefficients[2, 4])
		}
		
		# find the smallest p values
		index = sort(coeff, decreasing = FALSE, index.return = TRUE)$ix[1:5];
		
		# OLS
		# formula
		xname = paste0("V", index);
		formula = as.formula(paste("y ~ 1+", paste(xname, collapse = "+")));
		selected_data = as.data.frame(X[, index]);
		colnames(selected_data) = xname;
		pData = as.data.frame(XX[, index]);
		colnames(pData) = xname;
		
		fit = lm(formula, data = selected_data);
		conf_int = rbind(conf_int, confint(fit)[2:6,]);
		mse <- c(mse, mean((y - predict(fit, selected_data))^2));
		mspe <- c(mspe, mean((yy - predict(fit, pData))^2));
		setTxtProgressBar(pb, i);
	}
	list(mse = mse, mspe = mspe, conf_int = conf_int);
}