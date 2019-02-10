# generate data
mu = c(rep(3, 100), rep(0, 900)); # mean
Sigma = matrix(0, nrow = 1000, ncol = 1000); # covariance
Sigma[101:1000, 101:1000] = 0.1;
diag(Sigma) = 1;
library("MASS")
Z1 = mvrnorm(2000, mu, Sigma, empirical = TRUE);

Sigma = matrix(0, nrow = 1000, ncol = 1000);
diag(Sigma) = 1;
Z2 = mvrnorm(2000, mu, Sigma, empirical = TRUE);

fdr <- function(Z) {
  FDRs = rep(0, 2000);
  for(n in 1:2000) {
    # p values
    p = 2*pnorm(-abs(Z[n,]), 0, 1);
    
    # B-H with known pi0
    h = ncol(Z);
    q = 0.2;
    pi0 = 0.9;
    h0_hat = pi0*h;
    
    s = sort(p, index.return = TRUE);
    sorted_p = s$x;
    index = s$ix;
    
    i = 1;
    while(i <= h) {
      if(sorted_p[i] > i/h0_hat*q) {
        break;
      }
      i = i+1;
    }
    i = i-1;
    
    if(i == 0) {
      FDRs[n] = 0;
    }
    else {
      count = 0;
      for(j in 1:i) {
        if(index[j] > 100) # rejected null true
          count = count + 1;
      }
      FDRs[n] = count/i;
    }
  }
  return(FDRs);
}