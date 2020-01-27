require(robustbase)
# --------------------------------------------------------------------
# Author:
# *Enter your group number here, as well as names and student numbers*
# --------------------------------------------------------------------

## Use this code skeleton to implement the deterministic MCD and the plug-in
## robust regression estimator.  Please submit your implementation together
## with your report.

## IMPORTANT: Please do not change any function names and make sure that you
##            use the correct input and output as specified for each function.
##            This will simplify grading because each group's code will be
##            similarly structured.  However, feel free to add other arguments
##            to the function definitions as needed.



## Functions for initial estimators

# Input: the standardized data matrix z
# Output: the estimated covariance or correlation matrix
# Please do not use any other input or output for the initial estimators

# Fisher's transformation of the correlation coefficient
# correlation matrix based on the hyperbolic tangent
corHT <- function(z) {
  y = tanh(z) 
  correlation = cor(y, method = "pearson")
  return(correlation)
}

# spearman correlation matrix
corSpearman <- function(z) {
  n = dim(z)[1]
  p = dim(z)[2]
  
  #compute ranks of the columns
  ranks = matrix(NA, n, p)
  for (i in 1:p){
    ranks[,i] = rank(z[,i])
  }
  correlation = cor(ranks, method = "pearson")
  return(correlation)
}

# correlation matrix based on normal scores of the ranks
corNSR <- function(z) {
  n = dim(z)[1]
  p = dim(z)[2]
  #compute ranks of the columns
  ranks = matrix(NA, n, p)
  for (i in 1:p){
    ranks[ ,i] = rank(z[ ,i])
  }
  #compute normal scores from the ranks based on normal cumulative distribution
  normScores = qnorm((ranks - 1/3)/(n + 1/3))
  correlation = cor(normScores, method = "pearson")
  return(correlation)
}

# modified spatial sign covariance matrix
covMSS <- function(z) {
  n = dim(z)[1]
  p = dim(z)[2]
  norms = sqrt(rowSums(z^2)) #norms of rows
  k = as.matrix(z/matrix(rep(norms, p), n, p))
  covariance = 1/n*(t(k)%*%k)
  return(covariance)
}

# covariance matrix based on first step of BACON
covBACON1 <- function(z) {
  n = dim(z)[1]
  p = dim(z)[2]
  
  norms = sqrt(rowSums(z^2)) #norms of rows
  zSorted = z[order(norms),]
  c = ceiling(n/2)
  covariance = cov(zSorted[1:c, ], method = "pearson")
  return(covariance)
}

# raw OGK estimator of the covariance matrix with median and Qn
rawCovOGK <- function(z) {
  covariance = covOGK(z, sigmamu = s_mad)
  return(covariance$cov)
  # *enter your code here*
  # Hint: have a look at function covOGK() in package robustbase
}



## Main function for deterministic MCD algorithm

# Input:
# x ....... data matrix
# alpha ... proportion of observations to be used for subset size
# anything else you need

# Output
# A list with the following components:
# center ....... mean of the reweighted estimator
# cov .......... covariance matrix of the reweighted estimator
# weights ...... binary weights of the observations used for the reweighted
#                estimator (0 if outlier, 1 if used for estimation)
# raw.center ... mean of the raw estimator
# raw.cov ...... covariance matrix of the raw estimator
# best ......... indices of the observations in the best subset found by the
#                raw estimator
# any other output you want to return
#x <- Eredivisie28
covDetMCD <- function(x, alpha, ...) {
  #mahalo distance function
  n = dim(x)[1]
  p = dim(x)[2]
  
  mah_dist <- function(x, m_cov, v_mean){
    v_mah_dist <- vector()
    for(i in 1:nrow(x)){
      v_mah_dist[i] <- sqrt(as.matrix(x[i,] - v_mean) %*% inv(m_cov) %*% t(x[i,]-v_mean))
    }
    return(v_mah_dist)
  }
  
  qn_estimator <- function(x){
    constant = 2.21914
    h = floor(n/2)+1
    Qn = rep(NA, p)
    for (i in 1:p){
      variable = as.matrix(x[,i])
      rep_matrix = matrix(data = apply(variable, 2, function(x) rep(x, n)), ncol = ncol(variable)*n)
      distance = abs(mrx - t(mrx))
      distance = distance[upper.tri(distance, diag = FALSE)]
      distance = distance[order(distance)]
      
      Qn[i] = constant*distance[h*(h-1)/2];
    }
  }
  
  
  v_colmed <- colMedians(as.matrix(x))
  v_mean <- colMeans(m_z)
  qn? #have to define this
  m_z <- (x - v_colmed)/qn
  m_z <- x
  
  #initial estimates
  m_s <- list()
  m_s[[1]] <- corHT(m_z)
  m_s[[2]] <- corSpearman(m_z)
  m_s[[3]] <- corNSR(m_z)
  m_s[[4]] <- covMSS(m_z)
  m_s[[5]] <- covBACON1(m_z)
  m_s[[6]] <- rawCovOGK(m_z)
  
  #algorithm
  h <- as.integer(nrow(x)/2)
  while(h > 0){
    
    dist <- matrix()
    for(i in 1:6){
      m_e = eigen(m_s[[2]])
      m_b = Z %*% m_e
      m_l = # diagonal qn^2 of B
      m_s = m_e %*% m_l %*% t(m_e)
      v_estmu = chol2inv(m_s) %*% (colMedians(m_z %*% inv(chol2inv(m_s))))
      dist[,i] = mah_dist(x, m_s, v_estmu)
    }
    h <- h/2
    
  }
  # *enter your code here*
  #
  # Please note that the subset sizes for the MCD are not simply fractions of 
  # the number of observations in the data set, as discussed in the lectures.
  # You can use function h.alpha.n() from package robustbase to compute the 
  # subset size.
  
}



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

lmDetMCD <- function(x, y, alpha, ...) {
  # *enter your code here*
}

