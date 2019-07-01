
# Script to create observations:
# Inputs:
# n: number of observations per machine (int). Enter a string for mixed design (e.g "mixed")
# machine: enter "ideal" or "simple"
# model: enter 1 (int) for individual model and model 1, enter 2 (int) for model 2
# varXi: c_2^2, the subcases of the Ideal machines 2. enter 0.1 for 
#         subcase 3, 0.5 for subcase 2 and 0.9 for subcase 1. Set equal to 1 for Ideal machines 1.
# run: sets the seed for Gaussian fields generated. Does not work properly, so new Gaussian fields are enter
#     no matter what run is set to. Observations have the same locations for no matter what "run" is set to.

# Ouput:
# Dataframe:
# 100 first columns: locations of observations for 100 machines (x-values)
# column 101-200 : work of observations for 100 machines (z-values)
# pred.z: 100x5 matrix, where the rows represent the machine number 
#         and the columns represent the prediction points. This is the true work for the prediction points.
# n: number of observations per machine or a string if Mixed design.
# varXi: the value of varXi (input)

# Example of usage: 
# observations <- createObservations(60, "ideal", 1, 1, 1)
# this dataframe called observations, can now be used in the inference.R script (the function)
createObservations <- function(n, machine, model, varXi, run){
  
  ## inputs, variables that you can change ####
  N <- 100 # number of replicates
  obsvar0 <- 0.01^2  # variance of observation errors
  theta0 <- 0.65  # theta
  pred.x <- c(1, 2, 3, 6, 8)  # prediction points
  len.pred <-length(pred.x)
  
  ## creating observations for ideal machine ####
  if (identical(machine, "ideal")) {  # ideal machine
    
    ## SPDE parameters ####
    sigma0 <- 0.2/sqrt(0.5)
    range0 <- sqrt(2)
    alpha <- 2
    lambda0 <- alpha-1/2
    kappa0 <- sqrt(8*lambda0)/range0
    tau0 <- sqrt(gamma(lambda0)/(gamma(alpha)*sqrt(4*pi)*kappa0^3*sigma0^2))
    
    
    ## Create x and z, different or same number of observations for all machines ####
    if(is.character(n)){  # different number of observations for all machines
      n1 <- 60  # number of observations per replicates for machine 1 - 50
      n2 <- 5  # number of observations per replicates for machine 51 - 60
      n3 <- 5  # number of observations per replicates for machine 61 - 70
      n4 <- 5  # number of observations per replicates for machine 71 - 80
      n5 <- 6  # number of observations per replicates for machine 81 - 90
      n6 <- 5 # number of observations per replicates for machine 91 - 100
      
      for (i in 1:N){
        namex <- paste(sprintf("x%03d", i))
        namez <- paste(sprintf("z%03d", i))
        map <- if(i <= 50) 1 else ((i + ((10 - i) %% 10)) / 10 ) - 4  # maps i to 1,2,3,4,5,6 depending on how many observations the machine will have.
        switch(map,  # generates x and observation errors for the different machines
               {set.seed(i); temp.x <- sort(runif(n1,0.2,4)); temp.obs <- rnorm(n1,0,sd = sqrt(obsvar0))},
               {set.seed(i); temp.x <- sort(runif(n2,0.2,1.5)); temp.obs <- rnorm(n2,0,sd = sqrt(obsvar0))},
               {set.seed(i); temp.x <- sort(runif(n3,1.5,2.5)); temp.obs <- rnorm(n3,0,sd = sqrt(obsvar0))},
               {set.seed(i); temp.x <- sort(runif(n4,2.5,4)); temp.obs <- rnorm(n4,0,sd = sqrt(obsvar0))},
               {set.seed(i); temp.x <- sort(c(runif(n5/2,0.2,1), runif(n5/2,3,4))); temp.obs <- rnorm(n5,0,sd = sqrt(obsvar0))},
               {set.seed(i); temp.x <- sort(runif(n6,0.2,4)); temp.obs <- rnorm(n6,0,sd = sqrt(obsvar0))})
        assign(namex, temp.x)
        assign(namez, theta0*temp.x + temp.obs)  # z = theta*x + obs.error, discrepanc(y/ies) are added later 
      }
    }
    else{  # same number of observations for all machines
      for (i in 1:N) {
        namex <- paste(sprintf("x%03d", i))
        namez <- paste(sprintf("z%03d", i))
        set.seed(i)
        temp.x <- sort(runif(n,0.2,4))
        temp.obs <- rnorm(n,0,sd = sqrt(obsvar0))
        assign(namex, temp.x)
        assign(namez, theta0*temp.x + temp.obs)  # z = theta*x + obs.error, discrepanc(y/ies) are added later 
      }
    }

    
    ## Creating xi and adding to z ####
    pred.z <- matrix(rep(0, N*len.pred), nrow = N)
    for (i in 1:N){  # Gaussian process to create xi
      obs.x.pred.x <- sort(c(get(sprintf("x%03d", i)), pred.x))
      x.index.pred <- match(pred.x, obs.x.pred.x)
      D <- fields::rdist(obs.x.pred.x)
      sigma <- sigma0^2 * (1 + sqrt(3)*D / (range0/2) ) * exp(-sqrt(3)*D / (range0/2))
      
      namexi <- paste(sprintf("xi%03d", i))
      namez <- paste(sprintf("z%03d", i))
      
      set.seed(run*N+i)
      xi.pred <- sqrt(varXi)*Rfast::rmvnorm(1, mu = rep(0,times = length(get(sprintf("x%03d", i))) + len.pred), sigma = sigma)
      
      pred.z[i,] <- theta0*pred.x  + xi.pred[x.index.pred]  # getting z for predictions
      assign(namexi, xi.pred[-x.index.pred])  # getting xi of observations
      assign(namez, get(sprintf("z%03d", i)) + get(sprintf("xi%03d", i)))  # adding xi to z
    }
    ## Creating delta and adding to z if model 2 ####
    if (model == 2) {  
      #creating delta
      delta.x <- round(seq(0.2,8,0.01), 2)
      x.index.pred <- match(pred.x, delta.x)
      D <- fields::rdist(delta.x)
      sigma <- sigma0^2 * (1 + sqrt(3)*D / (range0/2) ) * exp(-sqrt(3)*D / (range0/2))
      set.seed(run)
      delta <- sqrt(1-varXi)*as.vector(Rfast::rmvnorm(1, mu = rep(0,times = length(delta.x)), sigma = sigma))
      slopes <- diff(delta)/diff(delta.x)
      
      pred.delta <- delta[x.index.pred]
      
      for (i in 1:100) {
        namedelta <- paste(sprintf("delta%03d", i))
        namez <- paste(sprintf("z%03d", i))
        temp.x <- get(sprintf("x%03d", i))
        
        map.x <- round(floor(temp.x*100)/100,2)
        x.index.map <- match(map.x, delta.x)
        map.z <- delta[x.index.map]
        pred.z[i,] <- pred.z[i,]  + pred.delta # getting z for predictions
        assign(namedelta, slopes[x.index.map]*(temp.x - map.x) + map.z)
        assign(namez, get(sprintf("z%03d", i)) + get(sprintf("delta%03d", i)))  # adding delta and xi to z
        }
    }
}
  
  
  ## creating observations for simple machine ####
  if (identical(machine, "simple")) {  # simple machine
    ## true process ####
    true_process <- function(x, theta, a)(theta*x)/(1+x/a)
    set.seed(run)
    a <- rnorm(N, 20, sd = 2)  # model discrepancy
        ## store true values for prediction points ####
    pred.z <- matrix(rep(0, N*len.pred), nrow = N)
    
    for (i in 1:N) {
      pred.z[i,] <- true_process(pred.x, theta0, a[i])
    }
    ## Create x and z, different or same number of observations for all machines ####
    if (is.character(n)) {  # different number of observations for all machines
      n1 <- 60  # number of observations per replicates for machine 1 - 50
      n2 <- 5  # number of observations per replicates for machine 51 - 60
      n3 <- 5  # number of observations per replicates for machine 61 - 70
      n4 <- 5  # number of observations per replicates for machine 71 - 80
      n5 <- 6  # number of observations per replicates for machine 81 - 90
      n6 <- 5 # number of observations per replicates for machine 91 - 100
      
      for (i in 1:N){
        namex <- paste(sprintf("x%03d", i))
        namez <- paste(sprintf("z%03d", i))
        map <- if(i <= 50) 1 else ((i + ((10 - i) %% 10)) / 10 ) - 4  # maps i to 1,2,3,4,5,6 depending on how many observations the machine will have.
        switch(map,  # generates x and observation errors for the different machines
               {set.seed(i); temp.x <- sort(runif(n1,0.2,4)); temp.obs <- rnorm(n1,0,sd = sqrt(obsvar0))},
               {set.seed(i); temp.x <- sort(runif(n2,0.2,1.5)); temp.obs <- rnorm(n2,0,sd = sqrt(obsvar0))},
               {set.seed(i); temp.x <- sort(runif(n3,1.5,2.5)); temp.obs <- rnorm(n3,0,sd = sqrt(obsvar0))},
               {set.seed(i); temp.x <- sort(runif(n4,2.5,4)); temp.obs <- rnorm(n4,0,sd = sqrt(obsvar0))},
               {set.seed(i); temp.x <- sort(c(runif(n5/2,0.2,1), runif(n5/2,3,4))); temp.obs <- rnorm(n5,0,sd = sqrt(obsvar0))},
               {set.seed(i); temp.x <- sort(runif(n6,0.2,4)); temp.obs <- rnorm(n6,0,sd = sqrt(obsvar0))})
        assign(namex, temp.x)
        assign(namez, true_process(temp.x, theta0, a[i]) + temp.obs)  # z = true process of simple machine + observation errror
      }
    } else {  # same number of observations for all machines
      for (i in 1:N) {
        namex <- paste(sprintf("x%03d", i))
        namez <- paste(sprintf("z%03d", i))
        set.seed(i)
        temp.x <- sort(runif(n,0.2,4))
        temp.obs <- rnorm(n,0,sd = sqrt(obsvar0))
        assign(namex, temp.x)
        assign(namez, true_process(temp.x, theta0, a[i]) + temp.obs)  # z = true process of simple machine + observation errror
      }
    }
  }
  
  
  ## return list with x, z and variance of xi####
  xnames <- znames <- c()
  observations <- vector(mode ="list", length = N*2 + 3)
  
  for (i in 1:N) {
    xnames <- c(xnames, sprintf("x%03d", i))
    znames <- c(znames, sprintf("z%03d", i))
    observations[[i]] <- get(sprintf("x%03d", i))
    observations[[N + i]] <- get(sprintf("z%03d", i))
  }
  observations[[2*N+1]] <- pred.z
  observations[[2*N+2]] <- n
  observations[[2*N+3]] <- varXi
  names(observations) <- c(xnames, znames, "pred.z","n","varXi")
  
  return(observations)
}
