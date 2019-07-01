
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

# prediction plots, data1 and data2 = dataframes from inference.R. 
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

# j for plotBox can only obtain value 1 or 2, for interpolation points and extrapolation points respectively.
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

printCoverageColSum12 <- function(data1, data2){
  N <- 100
  tableSingle <- table1 <- table2 <- matrix(rep(0, 5*N), nrow = N)
  for (i in 1:5){
    tableSingle[,i] <- data1[[sprintf("Coverage.single.%d", i)]][1:N]
    table1[,i] <- data1[[sprintf("Coverage.multiple.%d", i)]][1:N]
    table2[,i] <- data2[[sprintf("Coverage.multiple.%d", i)]][1:N]
  }
  ps <- sum(tableSingle)/500
  p1 <- sum(table1)/500
  p2 <- sum(table2)/500
  single <- c(colSums(tableSingle)/100, ps, ps-1.96*sqrt(ps*(1-ps)/500), ps+1.96*sqrt(ps*(1-ps)/500), pvalue(ps*500, 500))
  model1 <- c(colSums(table1)/100, p1, p1-1.96*sqrt(p1*(1-p1)/500), p1+1.96*sqrt(p1*(1-p1)/500), pvalue(p1*500, 500))
  model2 <- c(colSums(table2)/100, p2, p2-1.96*sqrt(p2*(1-p2)/500), p2+1.96*sqrt(p2*(1-p2)/500), pvalue(p2*500, 500))
  print("Coverage for individual model")
  print(round(100*single,1))
  print("Coverage for model 1")
  print(round(100*model1,1))
  print("Coverage for model 2")
  print(round(100*model2,1))
}

printCoverageVobs12 <- function(data1, data2, print1, print2){
  nDigits <- function(x) nchar( trunc( abs(x) ) )
  hspace <- function(x)if (nDigits(x) == 3) return(3) else return (5)
  N <- 100
  tableSingle <- table1 <- table2 <- matrix(rep(0, 5*N), nrow = N)
  for (i in 1:5){
    tableSingle[,i] <- data1[[sprintf("Coverage.single.%d", i)]]
    table1[,i] <- data1[[sprintf("Coverage.multiple.%d", i)]]
    table2[,i] <- data2[[sprintf("Coverage.multiple.%d", i)]]
  }
  
  single <- matrix(rep(0, 6*12), nrow = 6)
  model1 <- matrix(rep(0, 6*12), nrow = 6)
  model2 <- matrix(rep(0, 6*12), nrow = 6)
  single[1,1:5] <- colSums(tableSingle[1:50,])/50
  model1[1,1:5] <- colSums(table1[1:50,])/50
  model2[1,1:5] <- colSums(table2[1:50,])/50
  for (i in 2:6) {
    single[i,1:5] <- colSums(tableSingle[(31+i*10):(40+i*10),])/10
    model1[i,1:5] <- colSums(table1[(31+i*10):(40+i*10),])/10
    model2[i,1:5] <- colSums(table2[(31+i*10):(40+i*10),])/10
  }
  
  single[,6] <- rowSums(single[,1:5])/5
  model1[,6] <- rowSums(model1[,1:5])/5
  model2[,6] <- rowSums(model2[,1:5])/5
  mixed <- c(50,10,10,10,10,10)
  for (i in 1:5) {
    single[,i+6] <- pvalue(single[,i]*mixed,mixed)
    model1[,i+6] <- pvalue(model1[,i]*mixed,mixed)
    model2[,i+6] <- pvalue(model2[,i]*mixed,mixed)
  }
  single[,12] <- pvalue(single[,6]*mixed*5,mixed*5)
  model1[,12] <- pvalue(model1[,6]*mixed*5,mixed*5)
  model2[,12] <- pvalue(model2[,6]*mixed*5,mixed*5)
  # CI.single <- CI.model1 <- CI.model2 <- matrix(rep(0,times=6*2), nrow = 6)
  # CI.single[,1] <- single[,6] - 1.96*sqrt(single[,6]*(1-single[,6])/c(250,50,50,50,50,50))
  # CI.model1[,1] <- model1[,6] - 1.96*sqrt(model1[,6]*(1-model1[,6])/c(250,50,50,50,50,50))
  # CI.model2[,1] <- model2[,6] - 1.96*sqrt(model2[,6]*(1-model2[,6])/c(250,50,50,50,50,50))
  # CI.single[,2] <- single[,6] + 1.96*sqrt(single[,6]*(1-single[,6])/c(250,50,50,50,50,50))
  # CI.model1[,2] <- model1[,6] + 1.96*sqrt(model1[,6]*(1-model1[,6])/c(250,50,50,50,50,50))
  # CI.model2[,2] <- model2[,6] + 1.96*sqrt(model2[,6]*(1-model2[,6])/c(250,50,50,50,50,50))
  # single <- cbind(single, CI.single, pvalue(single[,6]*c(250,50,50,50,50,50),c(250,50,50,50,50,50)))
  # model1 <- cbind(model1, CI.model1, pvalue(model1[,6]*c(250,50,50,50,50,50),c(250,50,50,50,50,50)))
  # model2 <- cbind(model2, CI.model2, pvalue(model2[,6]*c(250,50,50,50,50,50),c(250,50,50,50,50,50)))
  print("Coverage table")
  for (i in 1:6) {
    cat(sprintf("%d & 
                %01d\\hspace{%dmm}(%01d) & %01d\\hspace{%dmm}(%01d) & 
                %01d\\hspace{%dmm}(%01d) & %01d\\hspace{%dmm}(%01d) &
                %01d\\hspace{%dmm}(%01d) & %01d\\hspace{%dmm}(%01d)", i,
                round(100*single[i,print1]), hspace(round(100*single[i,print1])), round(100*single[i,print1+6]),
                round(100*single[i,print2]), hspace(round(100*single[i,print2])), round(100*single[i,print2+6]),
                round(100*model1[i,print1]), hspace(round(100*model1[i,print1])), round(100*model1[i,print1+6]),
                round(100*model1[i,print2]), hspace(round(100*model1[i,print2])), round(100*model1[i,print2+6]),
                round(100*model2[i,print1]), hspace(round(100*model2[i,print1])), round(100*model2[i,print1+6]), 
                round(100*model2[i,print2]), hspace(round(100*model2[i,print2])), round(100*model2[i,print2+6])))
    cat("\\\\ \\hline \n")
                }
}

getCoverageIndividual <- function(data){
  N <- 100
  tableSingle <- matrix(rep(0, 5*N), nrow = N)
  for (i in 1:5){
    tableSingle[,i] <- data[[sprintf("Coverage.single.%d", i)]][1:N]
  }
  single <- c(colSums(tableSingle)/100, sum(tableSingle)/500)
  return(single)
}

getCoverageMultipleMd <- function(data){
  N <- 100
  tableMultiple <- matrix(rep(0, 5*N), nrow = N)
  for (i in 1:5){
    tableMultiple[,i] <- data[[sprintf("Coverage.multiple.%d", i)]][1:N]
  }
  
  model2 <- matrix(rep(0, 6*6), nrow = 6)
  model2[1,1:5] <- colSums(tableMultiple[1:50,])/50
  for (i in 2:6) {
    model2[i,1:5] <- colSums(tableMultiple[(31+i*10):(40+i*10),])/10
  }
  model2[,6] <- rowSums(model2[,1:5])/5
  model2 <- cbind(model2, pvalue(model2[,6]*c(250,50,50,50,50,50),c(250,50,50,50,50,50)))
  return(model2)
}