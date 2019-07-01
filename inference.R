# function to perform inference

# Input:
# observations: the dataframe produced with create_observations.R script
# model: enter 1 (int) for individual model and Model 1, enter 2 for Model 2

# Output: 
# Dataframe:
# It should somehow be obvious what it returns from the names of the list's entries.

# Example of usage:
# this is the way to apply the three models to the same set of observations
# data1 <- inference(observations, 1) # individual model and Model 1
# data2 <- inference(observations, 2) # Model 2
# Note that the "single" results returned if model is set to 2 (int), is not a part of the thesis.
# data1 and data2 can now be used in the plots.R script to obtain boxplots, prediction plots and coverage plots.

# It can be smart to save the data frames, and then it is possible to import the datasets later.
#write.table(data1, "enter_name1.csv", sep = ",", col.names = T, append = F, row.names = F)
#write.table(data2, "enter_name2.csv", sep = ",", col.names = T, append = F, row.names = F)

inference <- function(observations, model){
  library(INLA)  # INLA
  attach(observations)
  ## N, and prediction points ####
  N <- 100 # number of replicates
  pred.x <- c(1, 2, 3, 6, 8)  # prediction points
  len.pred <-length(pred.x)

  ## Mesh ####
  s <-c(seq(-16,0,0.2),seq(0.1,4.1, 0.01), seq(4.2,20, 0.2))
  mesh1 <- inla.mesh.1d(loc = s)
  
  ## SPDE ####
  if (model == 1) sigma0 <- 0.2/sqrt(0.5) else sigma0 <- 0.2
  range0 <- sqrt(2)
  alpha <- 2
  lambda0 <- alpha-1/2
  kappa0 <- sqrt(8*lambda0)/range0
  tau0 <- sqrt(gamma(lambda0)/(gamma(alpha)*sqrt(4*pi)*kappa0^3*sigma0^2))
  
  spde <- inla.spde2.matern(mesh=mesh1, alpha=2, constr = FALSE,
                            theta.prior.mean = c(log(tau0),log(kappa0)),
                            theta.prior.prec = c(0.1,0.1))
  
  ## Creating x, z and repl ####
  x <- c()
  z <- c()
  repl <- c()
  
  for (i in 1:N){
    x <- c(x, get(sprintf("x%03d", i)))
    z <- c(z, get(sprintf("z%03d", i)))
    repl <- c(repl, rep(i,length(get(paste(sprintf("x%03d", i))))))
  }
  
  ## Model 1, inference ####
  if (model == 1) {
    ## Prediction and estimation for multiple machines ####
    # making A, the projection matrix
    A.delta <- inla.spde.make.A(mesh1, loc = x, index = 1:length(x))
    A.xi <- inla.spde.make.A(mesh1, loc = x, index = 1:length(x), n.repl = N, repl = repl)
    A.pred.xi <- inla.spde.make.A(mesh=mesh1, loc=pred.x, index = rep(1:length(pred.x), times = N), 
                                  n.repl = N, repl = rep(1:N, each = length(pred.x)))
    
    # formula
    formula <- y ~ -1 + x + f(xi.field, model=spde, replicate=xi.field.repl)
    
    # mesh index
    mesh.index.xi <- inla.spde.make.index(name = "xi.field", n.spde = spde$n.spde, n.repl = N)
    
    # stacks
    stack.est <- inla.stack(data=list(y=z),
                            A=list(A.xi, 1),
                            effects=list(mesh.index.xi, data.frame(intercept = 1, x = x)),
                            tag="est")  #  estimation
    
    stack.pred.xi <- inla.stack(data=list(xi=NA),
                                A=list(A.pred.xi),
                                effects=list(mesh.index.xi),
                                tag="pred.xi")  # predict xi
    
    stack.pred.response <- inla.stack(data=list(y=NA),
                                      A=list(A.pred.xi, 1),
                                      effects=list(mesh.index.xi, data.frame(intercept = 1, x = rep(pred.x, times= N))),
                                      tag="pred.response")  # predict z
    
    join.stack <- inla.stack(stack.est, stack.pred.xi, stack.pred.response)
    
    # inla
    output.join <- inla(formula,
                        data=inla.stack.data(join.stack),
                        control.predictor=list(A=inla.stack.A(join.stack), compute=TRUE),
                        control.family = list(hyper = list(prec = list(prior = "loggamma", 
                                                                       param = c(181/19, 81/95000)))))
    
    # Indexing prediction results and storing results
    index.pred.xi <- inla.stack.index(join.stack, tag="pred.xi")$data
    index.pred.response <- inla.stack.index(join.stack, tag="pred.response")$data
    
    output.xi<-inla.spde2.result(inla=output.join, name = "xi.field", spde = spde, do.transform = TRUE)
    
    
    ## Prediction and estimation for individual machine ####
    ## Prediction and estimation for simple machine ####
    # A.pred
    A.pred.single <- inla.spde.make.A(mesh=mesh1, loc=pred.x)
    
    for (i in 1:N){
      # temp variables makes code more readable
      temp.x <- get(sprintf("x%03d", i))
      temp.z <- get(sprintf("z%03d", i))
      
      # A
      nameA <- paste(sprintf("A%03d", i))
      assign(nameA, inla.spde.make.A(mesh1, loc = temp.x))
      temp.A <- get(sprintf("A%03d", i))
      
      # formula
      x <- temp.x
      formula.single <- y ~ -1 + x + f(xi.field.single, model = spde)
      
      
      # mesh index
      namemesh.index.xi <- paste(sprintf("mesh.index.xi%03d", i))
      assign(namemesh.index.xi, inla.spde.make.index(name = "xi.field.single", n.spde = spde$n.spde))
      temp.mesh.index.xi <- get(sprintf("mesh.index.xi%03d", i))
      
      # stack
      namestack.est <- paste(sprintf("stack.est%03d", i))
      assign(namestack.est, inla.stack(data=list(y=temp.z),
                                       A=list(temp.A, 1),
                                       effects=list(temp.mesh.index.xi, data.frame(intercept = 1, x = temp.x)),
                                       tag=sprintf("est%03d", i))) #Estimation
      
      namestack.pred.xi <- paste(sprintf("stack.pred.xi%03d", i))
      assign(namestack.pred.xi, inla.stack(data=list(xi=NA), 
                                           A=list(A.pred.single),
                                           effects=list(temp.mesh.index.xi),
                                           tag=sprintf("pred.xi%03d", i))) #Predict xi
      
      namestack.pred.response <- paste(sprintf("stack.pred.response%03d", i))
      assign(namestack.pred.response, inla.stack(data=list(y=NA), 
                                                 A=list(A.pred.single, 1),
                                                 effects=list(temp.mesh.index.xi, data.frame(intercept = 1, x = pred.x)),
                                                 tag=sprintf("pred.response%03d", i))) #Predict response
      
      temp.stack.est <- get(sprintf("stack.est%03d", i))
      temp.stack.pred.xi <- get(sprintf("stack.pred.xi%03d", i))
      temp.stack.pred.response <- get(sprintf("stack.pred.response%03d", i))
      
      namejoin.stack <- paste(sprintf("join.stack%03d", i))
      assign(namejoin.stack, inla.stack(temp.stack.est, temp.stack.pred.xi, temp.stack.pred.response))
      temp.join.stack <- get(sprintf("join.stack%03d", i))
      
      # inla
      nameoutput.join <- paste(sprintf("output.join%03d", i))
      assign(nameoutput.join, inla(formula.single,
                                   data=inla.stack.data(temp.join.stack),
                                   control.predictor=list(A=inla.stack.A(temp.join.stack), compute=TRUE),
                                   control.family = list(hyper = list(prec = list(prior = "loggamma", 
                                                                                  param = c(181/19, 81/95000))))))
      
      
      # Storing results of SPDE parameters from estimation
      nameoutput.xi <- paste(sprintf("output.xi%03d", i))
      assign(nameoutput.xi, inla.spde2.result(inla=get(sprintf("output.join%03d", i)), name = "xi.field.single", spde = spde, do.transform = TRUE))
      
      # indexing latent and response prediction
      nameindex.pred.xi <- paste(sprintf("index.pred.xi%03d", i))
      nameindex.pred.response <- paste(sprintf("index.pred.response%03d", i))
      assign(nameindex.pred.xi, inla.stack.index(temp.join.stack, tag=sprintf("pred.xi%03d", i))$data)
      assign(nameindex.pred.response, inla.stack.index(temp.join.stack, tag=sprintf("pred.response%03d", i))$data)
    }
    
    
    
  } 
  ## Model 2, inference ####
  if (model == 2){
    ## Prediction and estimation for multiple machines ####
    # making A, the projection matrix
    A.delta <- inla.spde.make.A(mesh1, loc = x, index = 1:length(x))
    A.xi <- inla.spde.make.A(mesh1, loc = x, index = 1:length(x), n.repl = N, repl = repl)
    A.pred.delta <- inla.spde.make.A(mesh=mesh1, loc=pred.x, index = rep(1:length(pred.x), times = N))
    A.pred.xi <- inla.spde.make.A(mesh=mesh1, loc=pred.x, index = rep(1:length(pred.x), times = N), 
                                  n.repl = N, repl = rep(1:N, each = length(pred.x)))
    
    # formula
    formula <- y ~ -1 + x + f(delta.field, model = spde) + f(xi.field, model=spde, replicate=xi.field.repl)
    
    # mesh index
    mesh.index.xi <- inla.spde.make.index(name = "xi.field", n.spde = spde$n.spde, n.repl = N)
    mesh.index.delta <- inla.spde.make.index(name = "delta.field", n.spde = spde$n.spde)
    
    # stacks
    stack.est <- inla.stack(data=list(y=z),
                            A=list(A.delta, A.xi, 1),
                            effects=list(mesh.index.delta, mesh.index.xi, data.frame(intercept = 1, x = x)),
                            tag="est")  #  estimation
    
    stack.pred.delta <- inla.stack(data=list(delta=NA),
                                   A=list(A.pred.delta),
                                   effects=list(mesh.index.delta),
                                   tag="pred.delta")  # predict delta
    
    stack.pred.xi <- inla.stack(data=list(xi=NA),
                                A=list(A.pred.xi),
                                effects=list(mesh.index.xi),
                                tag="pred.xi")  # predict xi
    
    stack.pred.response <- inla.stack(data=list(y=NA),
                                      A=list(A.pred.delta, A.pred.xi, 1),
                                      effects=list(mesh.index.delta, mesh.index.xi, data.frame(intercept = 1, x = rep(pred.x, times= N))),
                                      tag="pred.response")  # predict z
    
    join.stack <- inla.stack(stack.est, stack.pred.delta, stack.pred.xi, stack.pred.response)
    
    # inla
    output.join <- inla(formula,
                        data=inla.stack.data(join.stack),
                        control.predictor=list(A=inla.stack.A(join.stack), compute=TRUE),
                        control.family = list(hyper = list(prec = list(prior = "loggamma", 
                                                                       param = c(181/19, 81/95000)))))
    
    # Indexing prediction results and storing results
    index.pred.delta <- inla.stack.index(join.stack, tag="pred.delta")$data
    index.pred.xi <- inla.stack.index(join.stack, tag="pred.xi")$data
    index.pred.response <- inla.stack.index(join.stack, tag="pred.response")$data
    
    output.delta<-inla.spde2.result(inla=output.join, name = "delta.field", spde = spde, do.transform = TRUE)
    output.xi<-inla.spde2.result(inla=output.join, name = "xi.field", spde = spde, do.transform = TRUE)
    
    
    ## Prediction and estimation for individual machine ####
    # A.pred
    A.pred.single <- inla.spde.make.A(mesh=mesh1, loc=pred.x)
    
    for (i in 1:N){
      # temp variables makes code more readable
      temp.x <- get(sprintf("x%03d", i))
      temp.z <- get(sprintf("z%03d", i))
      
      # A
      nameA <- paste(sprintf("A%03d", i))
      assign(nameA, inla.spde.make.A(mesh1, loc = temp.x))
      temp.A <- get(sprintf("A%03d", i))
      
      # formula
      x <- temp.x
      formula.single <- y ~ -1 + x + f(delta.field.single, model = spde) + f(xi.field.single, model = spde)
      
      
      # mesh index
      namemesh.index.delta <- paste(sprintf("mesh.index.delta%03d", i))
      assign(namemesh.index.delta, inla.spde.make.index(name = "delta.field.single", n.spde = spde$n.spde))
      namemesh.index.xi <- paste(sprintf("mesh.index.xi%03d", i))
      assign(namemesh.index.xi, inla.spde.make.index(name = "xi.field.single", n.spde = spde$n.spde))
      temp.mesh.index.delta <- get(sprintf("mesh.index.delta%03d", i))
      temp.mesh.index.xi <- get(sprintf("mesh.index.xi%03d", i))
      
      # stack
      namestack.est <- paste(sprintf("stack.est%03d", i))
      assign(namestack.est, inla.stack(data=list(y=temp.z),
                                       A=list(temp.A, temp.A, 1),
                                       effects=list(temp.mesh.index.delta, temp.mesh.index.xi, data.frame(intercept = 1, x = temp.x)),
                                       tag=sprintf("est%03d", i))) #Estimation
      
      namestack.pred.delta <- paste(sprintf("stack.pred.delta%03d", i))
      assign(namestack.pred.delta, inla.stack(data=list(delta=NA), 
                                              A=list(A.pred.single),
                                              effects=list(temp.mesh.index.delta),
                                              tag=sprintf("pred.delta%03d", i))) #Predict delta
      
      namestack.pred.xi <- paste(sprintf("stack.pred.xi%03d", i))
      assign(namestack.pred.xi, inla.stack(data=list(xi=NA), 
                                           A=list(A.pred.single),
                                           effects=list(temp.mesh.index.xi),
                                           tag=sprintf("pred.xi%03d", i))) #Predict xi
      
      namestack.pred.response <- paste(sprintf("stack.pred.response%03d", i))
      assign(namestack.pred.response, inla.stack(data=list(y=NA), 
                                                 A=list(A.pred.single, A.pred.single, 1),
                                                 effects=list(temp.mesh.index.delta, temp.mesh.index.xi, data.frame(intercept = 1, x = pred.x)),
                                                 tag=sprintf("pred.response%03d", i))) #Predict response
      
      temp.stack.est <- get(sprintf("stack.est%03d", i))
      temp.stack.pred.delta <- get(sprintf("stack.pred.delta%03d", i))
      temp.stack.pred.xi <- get(sprintf("stack.pred.xi%03d", i))
      temp.stack.pred.response <- get(sprintf("stack.pred.response%03d", i))
      
      namejoin.stack <- paste(sprintf("join.stack%03d", i))
      assign(namejoin.stack, inla.stack(temp.stack.est, temp.stack.pred.delta, temp.stack.pred.xi, temp.stack.pred.response))
      temp.join.stack <- get(sprintf("join.stack%03d", i))
      
      # inla
      nameoutput.join <- paste(sprintf("output.join%03d", i))
      assign(nameoutput.join, inla(formula.single,
                                   data=inla.stack.data(temp.join.stack),
                                   control.predictor=list(A=inla.stack.A(temp.join.stack), compute=TRUE),
                                   control.family = list(hyper = list(prec = list(prior = "loggamma", 
                                                                                  param = c(181/19, 81/95000))))))
      
      
      # Storing results of SPDE parameters from estimation
      nameoutput.delta <- paste(sprintf("output.delta%03d", i))
      nameoutput.xi <- paste(sprintf("output.xi%03d", i))
      assign(nameoutput.delta, inla.spde2.result(inla=get(sprintf("output.join%03d", i)), name = "delta.field.single", spde = spde, do.transform = TRUE))
      assign(nameoutput.xi, inla.spde2.result(inla=get(sprintf("output.join%03d", i)), name = "xi.field.single", spde = spde, do.transform = TRUE))
      
      # indexing latent and response prediction
      nameindex.pred.delta <- paste(sprintf("index.pred.delta%03d", i))
      nameindex.pred.xi <- paste(sprintf("index.pred.xi%03d", i))
      nameindex.pred.response <- paste(sprintf("index.pred.response%03d", i))
      assign(nameindex.pred.delta, inla.stack.index(temp.join.stack, tag=sprintf("pred.delta%03d", i))$data)
      assign(nameindex.pred.xi, inla.stack.index(temp.join.stack, tag=sprintf("pred.xi%03d", i))$data)
      assign(nameindex.pred.response, inla.stack.index(temp.join.stack, tag=sprintf("pred.response%03d", i))$data)
      
    }
  }
  ## Storing prediction results ####
  prediction.low.multiple <- prediction.mean.multiple <- prediction.high.multiple <- prediction.sd.multiple <- matrix(rep(0,N*len.pred), nrow = N)
  prediction.low.single <- prediction.mean.single <- prediction.high.single <- prediction.sd.single <- matrix(rep(0,N*len.pred), nrow = N)
  
  for (i in 1:N){
    for (j in 1:len.pred){
      prediction.low.multiple[i,j] <- output.join$summary.fitted.values[index.pred.response[(i-1)*len.pred + j], c("0.025quant")]
      prediction.mean.multiple[i,j] <- output.join$summary.fitted.values[index.pred.response[(i-1)*len.pred + j], c("mean")]
      prediction.high.multiple[i,j] <- output.join$summary.fitted.values[index.pred.response[(i-1)*len.pred + j], c("0.975quant")]
      prediction.sd.multiple[i,j] <- output.join$summary.fitted.values[index.pred.response[(i-1)*len.pred + j], c("sd")]
      prediction.low.single[i,j] <- get(sprintf("output.join%03d", i))$summary.fitted.values[get(sprintf("index.pred.response%03d",i))[j],c("0.025quant")]
      prediction.mean.single[i,j] <- get(sprintf("output.join%03d", i))$summary.fitted.values[get(sprintf("index.pred.response%03d",i))[j],c("mean")]
      prediction.high.single[i,j] <- get(sprintf("output.join%03d", i))$summary.fitted.values[get(sprintf("index.pred.response%03d",i))[j],c("0.975quant")]
      prediction.sd.single[i,j] <- get(sprintf("output.join%03d", i))$summary.fitted.values[get(sprintf("index.pred.response%03d",i))[j],c("sd")]
    }
  }
  
    # Counting number of points outside credible intervals (coverage percantage)
    coverage.multiple <- coverage.single <- matrix(rep(F,N*len.pred), nrow = N)
    for (i in 1:N){
      for (j in 1:len.pred){
        if (pred.z[i,j] <= prediction.high.multiple[i,j] & pred.z[i,j] >= prediction.low.multiple[i,j]){
          coverage.multiple[i,j] = T;
        }
        if (pred.z[i,j] <= prediction.high.single[i,j] & pred.z[i,j] >= prediction.low.single[i,j]){
          coverage.single[i,j] = T;
        }
      }
    }
    
    # Length of credible interval:
    CI.single <- prediction.high.single-prediction.low.single
    CI.multiple <- prediction.high.multiple-prediction.low.multiple
    
    ## MSE
    MSE.single <- (pred.z - prediction.mean.single)^2
    MSE.multiple <- (pred.z - prediction.mean.multiple)^2
    
    ## CRPS
    CRPS.z.single <- (pred.z - prediction.mean.single)/prediction.sd.single
    CRPS.z.multiple <- (pred.z - prediction.mean.multiple)/prediction.sd.multiple
    CRPS.single <- - prediction.sd.single*((1/sqrt(pi)) - 2*dnorm(CRPS.z.single) - CRPS.z.single*(2*pnorm(CRPS.z.single) - 1) )
    CRPS.multiple <- - prediction.sd.multiple*((1/sqrt(pi)) - 2*dnorm(CRPS.z.multiple) - CRPS.z.multiple*(2*pnorm(CRPS.z.multiple) - 1) )
  
  ## Storing parameter results ####
  # Storing parameters for multiple machines
  thetas <- output.join$summary.fixed$mean
  thetas0.025 <- output.join$summary.fixed$`0.025quant`
  thetas0.975 <- output.join$summary.fixed$`0.975quant`
  
  obsvars <- 1/output.join$summary.hyperpar[1,"mean"]
  obsvars0.025 <- 1/output.join$summary.hyperpar[1,"0.975quant"]
  obsvars0.975 <- 1/output.join$summary.hyperpar[1,"0.025quant"]
  
  if (model == 2){
    taus.delta <- exp(output.delta$summary.log.tau$mean)
    taus.delta0.025 <- exp(output.delta$summary.log.tau$`0.025quant`)
    taus.delta0.975 <- exp(output.delta$summary.log.tau$`0.975quant`)
    
    kappas.delta <- exp(output.delta$summary.log.kappa$mean)
    kappas.delta0.025 <- exp(output.delta$summary.log.kappa$`0.025quant`)
    kappas.delta0.975 <- exp(output.delta$summary.log.kappa$`0.975quant`)
  }
  
  taus.xi <- exp(output.xi$summary.log.tau$mean)
  taus.xi0.025 <- exp(output.xi$summary.log.tau$`0.025quant`)
  taus.xi0.975 <- exp(output.xi$summary.log.tau$`0.975quant`)
  
  kappas.xi <- exp(output.xi$summary.log.kappa$mean)
  kappas.xi0.025 <- exp(output.xi$summary.log.kappa$`0.025quant`)
  kappas.xi0.975 <- exp(output.xi$summary.log.kappa$`0.975quant`)
  
  # Storing parameters for individual machines
  thetas.single <- thetas0.025.single <- thetas0.975.single <- rep(0,N)
  obsvars.single <- obsvars0.025.single <- obsvars0.975.single <- rep(0,N)
  if(model == 2){
    taus.delta.single <- taus.delta0.025.single <-  taus.delta0.975.single <- rep(0,N)
    kappas.delta.single <- kappas.delta0.025.single <- kappas.delta0.975.single <- rep(0,N)
    
    for (i in 1:N){
      temp.output.delta <- get(sprintf("output.delta%03d",i))
      taus.delta.single[i] <- exp(temp.output.delta$summary.log.tau$mean)
      taus.delta0.025.single[i] <- exp(temp.output.delta$summary.log.tau$`0.025quant`)
      taus.delta0.975.single[i] <- exp(temp.output.delta$summary.log.tau$`0.975quant`)
      
      kappas.delta.single[i] <- exp(temp.output.delta$summary.log.kappa$mean)
      kappas.delta0.025.single[i] <- exp(temp.output.delta$summary.log.kappa$`0.025quant`)
      kappas.delta0.975.single[i] <- exp(temp.output.delta$summary.log.kappa$`0.975quant`)
    }
  }
  taus.xi.single <- taus.xi0.025.single <- taus.xi0.975.single <- rep(0,N)
  kappas.xi.single <- kappas.xi0.025.single <- kappas.xi0.975.single <- rep(0,N)
  
  for (i in 1:N){
    temp.output.join <- get(sprintf("output.join%03d",i))
    temp.output.xi <- get(sprintf("output.xi%03d",i))
    
    thetas.single[i] <- temp.output.join$summary.fixed$mean
    thetas0.025.single[i] <-  temp.output.join$summary.fixed$`0.025quant`
    thetas0.975.single[i] <-  temp.output.join$summary.fixed$`0.975quant`
    
    obsvars.single[i] <- 1/ temp.output.join$summary.hyperpar[1,"mean"]
    obsvars0.025.single[i] <- 1/ temp.output.join$summary.hyperpar[1,"0.975quant"]
    obsvars0.975.single[i] <- 1/ temp.output.join$summary.hyperpar[1,"0.025quant"]
    
    taus.xi.single[i] <- exp(temp.output.xi$summary.log.tau$mean)
    taus.xi0.025.single[i] <- exp(temp.output.xi$summary.log.tau$`0.025quant`)
    taus.xi0.975.single[i] <- exp(temp.output.xi$summary.log.tau$`0.975quant`)
    
    kappas.xi.single[i] <- exp(temp.output.xi$summary.log.kappa$mean)
    kappas.xi0.025.single[i] <- exp(temp.output.xi$summary.log.kappa$`0.025quant`)
    kappas.xi0.975.single[i] <- exp(temp.output.xi$summary.log.kappa$`0.975quant`)
  }
  ## Return ####
  if (model == 1){
      data_all <- data.frame("Coverage single" = coverage.single,
                             "Coverage multiple" = coverage.multiple,
                             "CI single" = CI.single,
                             "CI multiple" = CI.multiple,
                             "MSE single" = MSE.single,
                             "MSE multiple" = MSE.multiple,
                             "CRPS single" = CRPS.single,
                             "CRPS multiple" = CRPS.multiple,
                             "Observations" = pred.z,
                             "Prediction mean single" = prediction.mean.single,
                             "Prediction low single" = prediction.low.single,
                             "Prediction high single" = prediction.high.single,
                             "Prediction mean multiple" = prediction.mean.multiple,
                             "Prediction low multiple" = prediction.low.multiple,
                             "Prediction high multiple" = prediction.high.multiple,
                             "Theta multiple" = c(thetas, thetas0.025, thetas0.975, rep(NA, N-3)),
                             "Obsvar multiple" = c(obsvars, obsvars0.025, obsvars0.975, rep(NA, N-3)),
                             "Tau xi multiple" = c(taus.xi, taus.xi0.025, taus.xi0.975, rep(NA, N-3)),
                             "Kappa xi multiple" = c(kappas.xi, kappas.xi0.025, kappas.xi0.975, rep(NA, N-3)),
                             "Theta mean single" = thetas.single,
                             "Theta 0.025 single" = thetas0.025.single,
                             "Theta 0.975 single" = thetas0.975.single,
                             "Obsvar mean single" = obsvars.single,
                             "Obsvar 0.025 single" = obsvars0.025.single,
                             "Obsvar 0.975 single" = obsvars0.975.single,
                             "Tau xi mean single" = taus.xi.single,
                             "Tau xi 0.025 single" = taus.xi0.025.single,
                             "Tau xi 0.975 single" = taus.xi0.975.single,
                             "Kappa xi mean single" = kappas.xi.single,
                             "Kappa xi 0.025 single" = kappas.xi0.025.single,
                             "Kappa xi 0.975 single" = kappas.xi0.975.single)
  } else {
      data_all <- data.frame("Coverage single" = coverage.single,
                             "Coverage multiple" = coverage.multiple,
                             "CI single" = CI.single,
                             "CI multiple" = CI.multiple,
                             "MSE single" = MSE.single,
                             "MSE multiple" = MSE.multiple,
                             "CRPS single" = CRPS.single,
                             "CRPS multiple" = CRPS.multiple,
                             "Observations" = pred.z,
                             "Prediction mean single" = prediction.mean.single,
                             "Prediction low single" = prediction.low.single,
                             "Prediction high single" = prediction.high.single,
                             "Prediction mean multiple" = prediction.mean.multiple,
                             "Prediction low multiple" = prediction.low.multiple,
                             "Prediction high multiple" = prediction.high.multiple,
                             "Theta multiple" = c(thetas, thetas0.025, thetas0.975, rep(NA, N-3)),
                             "Obsvar multiple" = c(obsvars, obsvars0.025, obsvars0.975, rep(NA, N-3)),
                             "Tau delta multiple" = c(taus.delta, taus.delta0.025, taus.delta0.975, rep(NA, N-3)),
                             "Kappa delta multiple" = c(kappas.delta, kappas.delta0.025, kappas.delta0.975, rep(NA, N-3)),
                             "Tau xi multiple" = c(taus.xi, taus.xi0.025, taus.xi0.975, rep(NA, N-3)),
                             "Kappa xi multiple" = c(kappas.xi, kappas.xi0.025, kappas.xi0.975, rep(NA, N-3)),
                             "Theta mean single" = thetas.single,
                             "Theta 0.025 single" = thetas0.025.single,
                             "Theta 0.975 single" = thetas0.975.single,
                             "Obsvar mean single" = obsvars.single,
                             "Obsvar 0.025 single" = obsvars0.025.single,
                             "Obsvar 0.975 single" = obsvars0.975.single,
                             "Tau delta mean single" = taus.delta.single,
                             "Tau delta 0.025 single" = taus.delta0.025.single,
                             "Tau delta 0.975 single" = taus.delta0.975.single,
                             "Kappa delta mean single" = kappas.delta.single,
                             "Kappa delta 0.025 single" = kappas.delta0.025.single,
                             "Kappa delta 0.975 single" = kappas.delta0.975.single,
                             "Tau xi mean single" = taus.xi.single,
                             "Tau xi 0.025 single" = taus.xi0.025.single,
                             "Tau xi 0.975 single" = taus.xi0.975.single,
                             "Kappa xi mean single" = kappas.xi.single,
                             "Kappa xi 0.025 single" = kappas.xi0.025.single,
                             "Kappa xi 0.975 single" = kappas.xi0.975.single)
  }
  return(data_all)
}


