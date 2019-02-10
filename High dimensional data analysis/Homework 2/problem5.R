data <- readRDS("~/GitHub/Spring-19-courses/High dimensional data analysis/Homework 2/Singh2002.rds")
dx = data$X;
dy = data$y;

y = rep(0, length(dy));
for(i in 1:length(dy)) {
  if(dy[i] == "normal") {
    y[i] = 0;
  }
  else {
    y[i] = 1;
  }
}

x = matrix(rep(0, length(dx)), nrow = nrow(dx), ncol = ncol(dx));
for(i in 1:nrow(dx)) {
  for(j in 1:ncol(dx)) {
    x[i, j] = dx[i, j];
  }
}

