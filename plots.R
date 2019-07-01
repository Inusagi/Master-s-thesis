
# calculates pvalue, x = number of success, n =  number of trials.
pvalue <- function(x,n){
  ret <- rep(0,times=length(x))
  for (j in 1:length(x)) {
    if(x[j]/n[j] > 0.95){
      for(i in x[j]:n[j]){
        ret[j] <- ret[j] + dbinom(i,n[j],0.95)
      }
    } else if (x[j]/n[j] < 0.95){
      for(i in 0:x[j]){
        ret[j] <- ret[j] + dbinom(i,n[j],0.95)
      } 
    } else {
      ret[j] <- 1
    }
  }
  return(ret)
}

# prediction plots, 
# data1 and data2 = dataframes from inference.R, 
# scale: how much of y-values are seen. i or j = which prediction point (a number from 1-5)
plotPredictionCentralized12 <-function(data1, data2, scale, i){
  N <- 100
  len.pred <- 5
  pred.x <- c(1,2,3,6,8)
  
  prediction.mean.multiple <- data2[[sprintf("Prediction.mean.multiple.%d",i)]]
  ylim.Upper <- max((data2[[sprintf("Prediction.high.multiple.%d",i)]]-prediction.mean.multiple))*scale
  
  plot(1:N, data1[[sprintf("Prediction.mean.multiple.%d",i)]]-prediction.mean.multiple, type = "l", col = "limegreen", 
       ylim = c(-ylim.Upper, ylim.Upper),
       xlab = "Machine number", ylab = "Work - prediction mean of model 2")
  lines(1:N, data1[[sprintf("Prediction.low.multiple.%d",i)]]-prediction.mean.multiple, col = "limegreen", lty = 2)
  lines(1:N,data1[[sprintf("Prediction.high.multiple.%d",i)]]-prediction.mean.multiple, col = "limegreen", lty = 2)
  lines(1:N, data1[[sprintf("Prediction.mean.single.%d",i)]]-prediction.mean.multiple, col = "cornflowerblue")
  lines(1:N, data1[[sprintf("Prediction.low.single.%d",i)]]-prediction.mean.multiple, col = "cornflowerblue", lty = 2)
  lines(1:N, data1[[sprintf("Prediction.high.single.%d",i)]]-prediction.mean.multiple, col = "cornflowerblue", lty = 2)
  lines(1:N, data2[[sprintf("Prediction.mean.multiple.%d",i)]]-prediction.mean.multiple, col = "black")
  lines(1:N, data2[[sprintf("Prediction.low.multiple.%d",i)]]-prediction.mean.multiple, col = "black", lty = 2)
  lines(1:N,data2[[sprintf("Prediction.high.multiple.%d",i)]]-prediction.mean.multiple, col = "black", lty = 2)
  points(1:N, data1[[sprintf("Observations.%d",i)]]-prediction.mean.multiple, col = "red", pch = 4)
  
}

# coverage plots, 
# data1 and data2 = dataframes from inference.R, 
# scale: how much of y-values are seen. i or j = which prediction point (a number from 1-5)
plotPredictionCoverage12 <-function(data1, data2, scale, i){
  N <- 100
  len.pred <- 5
  pred.x <- c(1,2,3,6,8)
  prediction.mean.multiple <- data2[[sprintf("Prediction.mean.multiple.%d",i)]]
  ylim.Upper <- max((data2[[sprintf("Prediction.high.multiple.%d",i)]]-prediction.mean.multiple))*scale
  
  plot(1:N, data1[[sprintf("Prediction.mean.multiple.%d",i)]]-prediction.mean.multiple, type = "l", col = "white", 
       ylim = c(-ylim.Upper, ylim.Upper),
       xlab = "Machine number", ylab = "Work - prediction mean of model 2")
  lines(1:N, data1[[sprintf("Prediction.low.multiple.%d",i)]]-prediction.mean.multiple, col = "limegreen", lty = 1)
  lines(1:N,data1[[sprintf("Prediction.high.multiple.%d",i)]]-prediction.mean.multiple, col = "limegreen", lty = 1)
  lines(1:N, data1[[sprintf("Prediction.low.single.%d",i)]]-prediction.mean.multiple, col = "cornflowerblue", lty = 1)
  lines(1:N, data1[[sprintf("Prediction.high.single.%d",i)]]-prediction.mean.multiple, col = "cornflowerblue", lty = 1)
  lines(1:N, data2[[sprintf("Prediction.low.multiple.%d",i)]]-prediction.mean.multiple, col = "black", lty = 1)
  lines(1:N,data2[[sprintf("Prediction.high.multiple.%d",i)]]-prediction.mean.multiple, col = "black", lty = 1)
  points(1:N, data1[[sprintf("Observations.%d",i)]]-prediction.mean.multiple, col = "red", pch = 4)
}

# boxplots, 
# data1 and data2 = dataframes from inference.R, 
# scale: how much of y-values are seen. 
# j for plotBox can only obtain value 1 or 2, for interpolation points and extrapolation points respectively.
# accuracyMeasure "CRPS" or "MSE" or "CI"
plotBox12 <- function(data1, data2, accuracyMeasure, scale, j){
  len.pred <- 5
  accuracy.measure.single <- accuracy.measure.1 <- accuracy.measure.2 <-
    matrix(0,length(data1[[sprintf("%s.single.%d",accuracyMeasure,1)]])*len.pred, ncol = len.pred)
  
  ylim1 <- ylim2 <- 0
  for (i in 1:5){
    accuracy.measure.single[,i] <- data1[[sprintf("%s.single.%d", accuracyMeasure, i)]]
    accuracy.measure.1[,i] <- data1[[sprintf("%s.multiple.%d", accuracyMeasure, i)]]
    accuracy.measure.2[,i] <- data2[[sprintf("%s.multiple.%d", accuracyMeasure, i)]]
  }
  
  yBound1 <- max(accuracy.measure.single[,1:3],accuracy.measure.1[,1:3],accuracy.measure.2[,1:3])*scale
  yBound2 <- max(accuracy.measure.single[,4:5],accuracy.measure.1[,4:5],accuracy.measure.2[,4:5])*scale
  
  if (j == 1) {
    boxplot(accuracy.measure.single[,1], accuracy.measure.1[,1], accuracy.measure.2[,1],
            accuracy.measure.single[,2], accuracy.measure.1[,2], accuracy.measure.2[,2],
            accuracy.measure.single[,3], accuracy.measure.1[,3], accuracy.measure.2[,3], 
            border = rep(c("cornflowerblue", "limegreen", "black"), times = 3), ylim = c(0, yBound1),
            names = rep(seq(1,3,1),each = 3), xlab = "x-values of interpolation points")
  } else {
    boxplot(accuracy.measure.single[,4], accuracy.measure.1[,4], accuracy.measure.2[,4], 
            accuracy.measure.single[,5], accuracy.measure.1[,5], accuracy.measure.2[,5], 
            border = rep(c("cornflowerblue", "limegreen", "black"), times = 2 ), names = rep(c(6, 8), each = 3),
            xlab = "x-values of extrapolation points", ylim = c(0, yBound2))
  }
}

# prediction plots, 
# data1 and data2 = dataframes from inference.R, 
# scale: how much of y-values are seen. i or j = which prediction point (a number from 1-5)
# accuracyMeasure "CRPS" or "MSE" or "CI"
plotBoxVObs12 <- function(data1, data2, accuracyMeasure, scale, j){
  len.pred <- 5
  pred.x <- c(1,2,3,6,8)
  accuracy.measure.single <- accuracy.measure.multiple1 <- accuracy.measure.multiple2 <- 
    matrix(0,length(data1[[sprintf("%s.single.%d",accuracyMeasure,1)]])*len.pred, ncol = len.pred)
  
  for (i in 1:5){
    accuracy.measure.single[,i] <- data1[[sprintf("%s.single.%d", accuracyMeasure, i)]]
    accuracy.measure.multiple1[,i] <- data1[[sprintf("%s.multiple.%d", accuracyMeasure, i)]]
    accuracy.measure.multiple2[,i] <- data2[[sprintf("%s.multiple.%d", accuracyMeasure, i)]]
  }
  
  
  yBound <- max(accuracy.measure.single[,j],accuracy.measure.multiple1[,j],accuracy.measure.multiple2[,j])*scale
  boxplot(accuracy.measure.single[1:50,j], accuracy.measure.multiple1[1:50,j], accuracy.measure.multiple2[1:50,j],
          accuracy.measure.single[51:60,j], accuracy.measure.multiple1[51:60,j], accuracy.measure.multiple2[51:60,j], 
          accuracy.measure.single[61:70,j], accuracy.measure.multiple1[61:70,j], accuracy.measure.multiple2[61:70,j], 
          accuracy.measure.single[71:80,j], accuracy.measure.multiple1[71:80,j], accuracy.measure.multiple2[71:80,j],
          accuracy.measure.single[81:90,j], accuracy.measure.multiple1[81:90,j], accuracy.measure.multiple2[81:90,j],
          accuracy.measure.single[91:100,j], accuracy.measure.multiple1[91:100,j], accuracy.measure.multiple2[91:100,j],
          border = rep(c("cornflowerblue", "limegreen", "black"), times = 6), names = floor(seq(1,6.7,1/3)), 
          xlab = "Case", ylim = c(0, yBound))
  
}
