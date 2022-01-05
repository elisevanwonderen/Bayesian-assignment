
####################
#  Model results   #
####################

# Function to compute means, sds, credible intervals and the MC error:
get_estimates <- function(x) c(mean=mean(x), sd = sd(x),    # means and sds
                               MC_error = sd(x)/sqrt(length(x)),   # MC error
                               quant=quantile(x, probs = c(.025, .975))) # credible intervals

##################
#      DIC       #
##################

# Function to calculate the DIC. For "samples", a dataframe or matrix should be given as input
# that only contains the parameter estimates.
getDIC <- function(samples, y, x1, x2){

  # For Dhat, we obtain the mean parameter values, and then we sum the -2*loglikelihoods of
  # each y-value given the parameter values:
  theta.bar <- apply(samples, 2, mean) # obtain the mean parameter values
  Dhat <- -2*sum(dnorm(y, mean = theta.bar["b0"] + theta.bar["b1"]*x1 + theta.bar["b2"]*x2,
                       sd = sqrt(theta.bar["var"]), log = T))

  # For Dbar, instead of using the mean parameter values, we calculate a separate -2*loglikelihood
  # per sample, and take the mean of those -2*loglikelihoods:
  logliks <- rep(NA, nrow(samples)) # storage space for the log-likelihoods per sample
  for(i in 1:nrow(samples)){
    logliks[i] <- -2*sum(dnorm(y, mean = samples[i, "b0"] + samples[i, "b1"]*x1 + samples[i, "b2"]*x2,
                               sd = sqrt(samples[i, "var"]), log = T))
  }
  Dbar <- mean(logliks)

  # Now, we can calculate the DIC:
  DIC <- Dhat + 2*(Dbar - Dhat)
  return(DIC)
}
