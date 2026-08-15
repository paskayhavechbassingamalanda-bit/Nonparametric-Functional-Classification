###############################################################
# Data Generation Process
#
# Master's Thesis:
# Nonparametric Classification for Spatial Functional Data
# under Random Censoring
#
# Author: Paska Bassinga Malanda
###############################################################
#------------------------------------------------------------
# Function:
# Simulates spatial functional data under random censoring
#
# Arguments
# n1, n2 : spatial grid dimensions
# a       : spatial dependence parameter
# sigma   : standard deviation of the spatial noise
#
# Returns
# Matrix containing:
# X      : functional observations
# Y      : response
# delta  : censoring indicator
# Y_G    : transformed response
#------------------------------------------------------------
library(fda) 
library(fBasics)
library(fda.usc)
library(OmicsPLS)
library(fields)
library(geoR)
library(MultiRNG)
library(Matrix)


data2 = function(n1, n2, a, sigma) {
  
  nxn <- n1 * n2
  sigma2 <- sigma^2
  
  argvals = seq(0, 1, 1/364) ## Functional grid
  p = length(argvals)
  nbasis <- 19
  
  ## Spatial grid
  step <- 1
  x1 <- seq(1, n1, step)
  x2 <- seq(1, n2, step)
  coordinates <- expand.grid(x1, x2) ## Spatial grid
  
  mean.coef <- rep(0, nxn)
  distance <- dist(coordinates)
  distsite <- as.matrix(distance)
  
  Dep <- t(matrix(0, nxn, nxn))
  Dep <- exp(-a * (distsite))
  Dep <- apply(Dep, 2, mean)
  
  ###--------------------------------------------------------------------------------
  ###------------------------FUNCTIONAL BASIS----------------------------------------
  ###--------------------------------------------------------------------------------
  
  rangeval <- range(argvals)
  bf.basis <- create.fourier.basis(
    rangeval,
    nbasis,
    period = diff(rangeval)
  )
  
  Phi <- eval.basis(argvals, bf.basis, Lfdobj = 0)
  
  ###--------------------------------------------------------------------------------
  ###-------------------GENERATION OF FUNCTIONAL DATA-------------------------------
  ###--------------------------------------------------------------------------------
  
  X <- array(0, dim = c(nxn, p))
  
  coef_array <- matrix(
    0,
    nrow = nxn,
    ncol = nbasis
  )
  
  for (k in 1:nbasis) {
    
    covariance.coef <- cov.spatial(
      distsite,
      cov.model = "exponential",
      cov.pars = c(0.5^(0.5 * (k - 1)), 1/a)
    )
    
    covariance.coef <- as.matrix(covariance.coef) +
      diag(0.5^(k - 1), nxn)
    
    coef_vect <- mvrnorm(
      1,
      mu = mean.coef,
      covariance.coef
    )
    
    coef_array[, k] <- coef_vect
  }
  
  x <- coef_array %*% t(Phi)
  X <- Dep * x
  X = matrix(X, nxn, p)
  
  
  ###--------------------------------------------------------------------------------
  ###--------------------------------EPSILON-----------------------------------------
  ###--------------------------------------------------------------------------------
  
  covariance.coef1 <- cov.spatial(
    distance,
    cov.model = "exp",
    cov.pars = c(sigma, 1/a)
  )
  
  covariance.coef1 <- as.matrix(covariance.coef1) +
    diag(sigma2, nxn, nxn)
  
  covariance.coef1 = makePositiveDefinite(covariance.coef1)
  
  epsilon = mvrnorm(
    1,
    mean.coef,
    covariance.coef1
  )
  
  
  ###--------------------------------------------------------------------------------
  ###--------------------------RESPONSE VARIABLE-------------------------------------
  ###--------------------------------------------------------------------------------
  
  s = numeric(nxn)
  
  for (j in 1:p) {
    s = s + X[, j] * X[, j]
  }
  
  Y = pi * s / p + epsilon
  
  
  C <- rexp(
    nxn,
    rate = 0.2
  ) 
  # Censoring time 
  # (0.1 = 18%, 0.34 = 36%, 0.50 = 45%, 0.533 = 50%,
  #  0.60 = 53%, 0.90 = 63%, 1.5 = 72%, 2.5 = 81%,
  #  3 = 83%, 3.60 = 85%)
  
  C <- Dep * C
  
  Z <- pmin(Y, C)  # Observed time
  
  delta <- as.numeric(Y <= C)  
  # Censoring indicator (1 = uncensored, 0 = censored)
  
  
  # Calculate the quantile for the truncation point t0
  t0 <- quantile(Z, 0.95)
  
  # Create the empirical cumulative distribution function
  F_empirique <- ecdf(C)
  
  
  # Define the integrand
  integrale <- function(s) {
    1 / (1 - F_empirique(s))
  }
  
  
  # Numerical integration using the trapezoidal rule
  integrate_trapezoid <- function(func, lower, upper, n = 1000) {
    
    x <- seq(
      lower,
      upper,
      length.out = n
    )
    
    y <- func(x)
    
    h <- (upper - lower) / (n - 1)
    
    integral <- (h / 2) * (
      y[1] +
        2 * sum(y[2:(n - 1)]) +
        y[n]
    )
    
    return(integral)
  }
  
  
  # Apply the transformation only to censored observations
  Y_G <- numeric(nxn)
  
  for (i in 1:(nxn)) {
    
    top <- min(Z[i], t0)
    
    if (delta[i] == 0) {  
      # If the observation is censored
      
      Y_G[i] <- integrate_trapezoid(
        integrale,
        lower = 0,
        upper = top,
        n = 1000
      )
      
    } else {  
      # If the observation is not censored
      
      Y_G[i] <- Y[i]
    }
  }
  
  
  # Check the censoring percentage
  censorship_rate <- mean(delta == 0) * 100
  
  print(
    paste(
      "Censoring percentage:",
      censorship_rate,
      "%"
    )
  )
  
  
  Dn = cbind(
    X,
    Y,
    delta,
    Y_G
  )
  
  return(Dn)
  
}





