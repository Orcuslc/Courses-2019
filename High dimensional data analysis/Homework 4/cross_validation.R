source("lasso.R");

lasso.cv <- function(X, y, lambdas, nfolds = 10) {
	# normalize X
	#normalize <- function(x) return((x-mean(x))/sd(x));
	#X = apply(X, 2, normalize);
	
	# lambdas
	#lambdas = 10^seq(-2, 0, length = 100);
	
	# divide the data into nfolds 
	fold = floor(nrow(X)/nfolds);
	index = c((2:nfolds-1)*fold, nrow(X));
	
	# init cv error
	CV_error = rep(0, length(lambdas));
	
	# compute solution path
	for(k in 1:length(lambdas)) {
		for(i in 1:nfolds) {
			if(i == 1) {
				train_index = (index[1]+1):index[nfolds];
				test_index = 1:index[1];
			}
			else if(i == nfolds) {
				train_index = 1:index[nfolds-1];
				test_index = (index[nfolds-1]+1):index[nfolds];
			}
			else {
				train_index = c(1:index[i-1], (index[i]+1):index[nfolds]);
				test_index = (index[i-1]+1):index[i];
			}
			X_train = X[train_index, ];
			X_test = X[test_index, ];
			y_train = y[train_index];
			y_test = y[test_index];
			
			beta = lasso.coordinate_descent(X_train, y_train, lambdas[k]);
			#beta = glmnet(X_train, y_train, lambda = lambdas[k])$beta;
			y_pred = X_test %*% beta;
			pred_error = mean((y_test - y_pred)^2);

			CV_error[k] = CV_error[k] + pred_error;
		}
	}
	
	CV_error = CV_error/nfolds;
	return(lambdas[which.min(CV_error)]);
}