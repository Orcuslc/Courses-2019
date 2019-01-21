source("exercise1.2.R");

res_p5 = ex1.2(p = 5);
res_p50 = ex1.2(p = 50);
res_p500 = ex1.2(p = 500);
res_p5000 = ex1.2(p = 5000);

p = c(5, 50, 500, 5000);
mse = c(mean(res_p5$mse), mean(res_p50$mse), mean(res_p500$mse), mean(res_p5000$mse));
mspe = c(mean(res_p5$mspe), mean(res_p50$mspe), mean(res_p500$mspe), mean(res_p5000$mspe));
plot(p, mse, xlab = "p", ylab = "", type = "b", col = "red", log="x", xaxt = "n", yaxt = "n");
axis(side = 1, at = p, labels = p);
par(new=TRUE);
plot(p, mspe, xlab = "p", ylab = "", type = "b", col = "blue", log = "x", xaxt = "n", yaxt = "n", ylim = c(0.3, 1.4), lty = 2);
axis(side = 2, at = c(0.3, 0.5, 0.7, 0.9, 1.1, 1.3), labels = c(0.3, 0.5, 0.7, 0.9, 1.1, 1.3));
legend(5, 0.5, legend = c("MSE", "MSPE"), col = c("red", "blue"), lty = 1:2);
title("MSE/MSPE vs. p");

conf_int = rbind(apply(res_p5$conf_int, 2, mean), apply(res_p50$conf_int, 2, mean), apply(res_p500$conf_int, 2, mean), apply(res_p5000$conf_int, 2, mean))
plot.new();
require(plotrix);
plotCI(c(1, 2, 3, 4), apply(conf_int, 1, mean), ui = conf_int[, 2], li = conf_int[, 1], xaxt = "n", xlab = "p", ylab = "confidence interval");
axis(side = 1, at = c(1, 2, 3, 4), labels = c(5, 50, 500, 5000));
title("confidence interval vs. p")
