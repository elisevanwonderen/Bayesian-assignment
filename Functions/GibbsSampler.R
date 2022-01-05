############################################
#        GIBBS SAMPLER WITH M-H STEP       #
############################################

# Gibbs sampler for a regression analysis with 2 predictors, and a M-H step for b2.
# Normal priors are used for parameters b0-b2, and an inverse gamma prior for the variance.
# The user can specify prior means (mu) and variances (tau) for parameters b0-b2, and prior
# hyperparameters alpha and beta for the variance. The default values represent uninformative priors.
# We also need values for the mean (m.prop) and SD (sd.prop) for the proposal distribution for b2.
# The user can also specify the number of chains and different starting values for each chain.
# In addition, the default number of iterations is set to 10.000 with a burn-in period of 2000.
myGibbsSampler <- function(y, x1, x2, mu00 = 0, mu10 = 0, mu20 = 0,
                           tau00 = 1000, tau10 = 1000, tau20 = 1000,
                           alpha = 0.001, beta =0.001, m.prop, sd.prop,
                           inits, n.chains = 2, n.iter = 10000, burnin = 2000){

  set.seed(123)

  # Function to calculate the log probability of a value from the conditional
  # posterior density for the M-H step. We use the log to prevent values becoming infinite.
  log.p.post <- function(b2, b0, b1, var, mu20, tau20, x1, x2){
    mean <- ((sum(x2*(y-b0-b1*x1))/var) + (mu20/tau20)) / ((sum(x2^2)/var) + (1/tau20))
    sd <- sqrt(1/((sum(x2^2)/var) + (1/tau20)))
    log.prob <- dnorm(b2, mean, sd, log = T)
    return(log.prob)
  }

  # create storage space for parameter values (for all chains)
  samples <- as.data.frame(matrix(nrow = n.chains*n.iter, ncol = 9))
  colnames(samples) <- c("b0", "b1", "b2", "var", "accept", "log.post", "log.prop", "R", "u")
  seq <- 1:n.iter  # rows to be filled for the first chain

  n <- length(y) # sample size

  for(k in 1:n.chains){ # repeat the sequence below for each chain

    # extract starting values:
    # the first value for b0 is computed based on the other parameters
    # so we don't need a starting value for b0
    b1 <- as.numeric(inits[[k]]['b1'])
    b2 <- as.numeric(inits[[k]]['b2'])
    var <- as.numeric(inits[[k]]['var'])

    # Calculate the parameter values in each iteration as a function of
    # the prior parameters, the data, and the current values of the other
    # model parameters
    for(i in seq){
      b0 <- rnorm(1, mean = ((sum(y-b1*x1-b2*x2)/var) + (mu00/tau00)) / ((n/var) + (1/tau00)),
                  sd = sqrt(1/ ((n/var) + (1/tau00))))
      b1 <- rnorm(1, mean = ((sum(x1*(y-b0-b2*x2))/var) + (mu10/tau10)) / ((sum(x1^2)/var) + (1/tau10)),
                  sd = sqrt(1/((sum(x1^2)/var) + (1/tau10))))
      b.star <- rnorm(1, m.prop, sd.prop)  # sample value from proposal density
      log.prop <- dnorm(b2, m.prop, sd.prop, log = T) - dnorm(b.star, m.prop, sd.prop, log = T)
      log.post <-  log.p.post(b.star, b0, b1, var, mu20, tau20, x1, x2) -
        log.p.post(b2, b0, b1, var, mu20, tau20, x1, x2)
      R <- log.post + log.prop  # log probability of acceptance
      u <- log(runif(1, min = 0, max = 1))  # sample value from uniform distribution and take the log
      if(R >= u){
        b2 <- b.star    # accept proposed value if R is greater than or equal to u
        accept <- 1
      } else {
        b2 <- b2       # if not, retain old value
        accept <- 0
      }
      var <- 1/(rgamma(1, shape = (n/2) + alpha, rate = (sum((y-(b0+b1*x1+b2*x2))^2)/2) + beta))
      samples[i,] <- c(b0, b1, b2, var, accept, log.post, log.prop, R, u)  # save values
    }

    seq <- (k*n.iter+1):((k+1)*n.iter)  # rows to be filled for the next chain
  }
  samples$chain <- rep(1:n.chains, each = n.iter) # specify which chain the samples belong to
  samples$iteration <- rep(1:n.iter, times = n.chains) # add variable indicating the iteration
  samples <- samples[!(samples$iteration %in% c(1:burnin)),] # remove burn-in iterations
  return(samples)
}
