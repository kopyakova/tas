## Function for regression based on the deterministic MCD

# Input:
# x ........ matrix of explanatory variables
# y ........ response variable
# alpha .... proportion of observations to be used for the subset size in the 
#            MCD estimator
# anything else you need

# Output
# A list with the following components:
# coefficients .... estimated regression coefficients based on the reweighted
#                   deterministic MCD covariance matrix
# fitted.values ... fitted values for all observations in the data
# residuals ....... residuals for all observations in the data
# MCD ............. R object for the deterministic MCD (entire output from
#                   function covDetMCD())
# any other output you want to return

lmDetMCD <- function(x, y, alpha = 0.5, ...) {
  n = length(x)
  data <- cbind(x, y)
  
  #get MCD estimates
  MCD = covMcd(data, alpha, nsamp = "deterministic") 
  mu_x = MCD$center[1]
  mu_y = MCD$center[2]
  sigma_yy = MCD$cov[2,2]
  sigma_xx = MCD$cov[1,1]
  sigma_xy = MCD$cov[1,2]
  
  #calculate coefficients
  beta = solve(sigma_xx) %*% sigma_xy
  intercept = mu_y - mu_x*beta
  coefficients = c(intercept, beta)
  
  #calculate predicted / fitted values and residuals
  fitted.values = cbind(rep(1, n),x) %*% coefficients
  residuals = y - fitted.values
  
  return(list(coefficients, fitted.values, residuals, MCD))
}