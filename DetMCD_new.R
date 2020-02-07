library(robustbase)
library(expm)
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
  correlation = cor(z, method = "spearman")
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
  index = norms > .Machine$double.eps
  k = z
  k[index,] = z[index,] / norms[index]
  covariance = t(k) %*% k
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
  covariance = covOGK(z, sigmamu = s_Qn)
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

covDetMCD <- function(x, alpha = 0.5, maxiter = 100, delta = 0.025) {
  
  n = dim(x)[1]
  p = dim(x)[2]

  #Qn scale estimator function
  qn_estimator <- function(x){
    constant = 2.21914
    n = dim(x)[1]
    p = dim(x)[2]
    h = floor(n/2)+1
    Qn = rep(NA, p)
    for (i in 1:2){
      variable = as.matrix(x[,i])
      rep_matrix = matrix(data = apply(variable, 2, function(x) rep(x, n)), ncol = ncol(variable)*n)
      distance = abs(rep_matrix - t(rep_matrix))
      distance = distance[upper.tri(distance, diag = FALSE)]
      distance = distance[order(distance)]
      
      Qn[i] = constant*distance[h*(h-1)/2];
    }
    return(Qn)
  }
 
  #get medians of the columns
  v_colmed <-  apply(x, 2L, median)
  #get scale qn
  qn <- qn_estimator(x)
  #standardize the data
  m_z <- cbind(
    (x[,1]-v_colmed[1])*(1/qn[1]),
    (x[,2]-v_colmed[2])*(1/qn[2]))

  
  #initial estimates
  m_s <- list()

  m_s[[1]] <- corHT(m_z)
  m_s[[2]] <- corSpearman(m_z)
  m_s[[3]] <- corNSR(m_z)
  m_s[[4]] <- covMSS(m_z) 
  m_s[[5]] <- covBACON1(m_z)
  m_s[[6]] <- rawCovOGK(m_z) 

  #initialize data structures
  h <- h.alpha.n(alpha,n,p)
  dist  <- matrix(data = 0, ncol = 6, nrow = n)
  dist2 <- matrix(data = 0, ncol = 6, nrow = h)
  
  detnew <- matrix(data=0, ncol =6, nrow = maxiter)
  estmu2 <- list()
  m_s2 <- list()
  lowdet <- vector()

  for(i in 1:6){
    
    #estimate initial distance according to Hubert et al. (2012)
    m_e = eigen(m_s[[i]],  symmetric=TRUE)$vectors
    m_b = m_z %*% m_e
    m_l = diag(qn_estimator(m_b)^2) #qn^2 doesn't make a differen
    m_s2[[i]] = m_e %*% m_l %*% t(m_e)
    m_sqrt = sqrtm(m_s2[[i]])
    temp = m_z %*% solve(m_sqrt)
    v_estmu = m_sqrt %*% apply(temp, 2L, median)
    
    dist[,i] =  sqrt(mahalanobis(x = m_z, cov =  m_s2[[i]], center = v_estmu))
  
    #c - steps
    det_old <- Inf
     for (k in 1:maxiter){
       
      #select sample
      selected_sample <- x[sort.list(dist[,i])[1:h],]
      
      #get new estimates
      estmu2[[i]] <- apply(selected_sample, 2, mean)
      v_estmu     <- estmu2[[i]] 
      m_s2[[i]]   <- cov(selected_sample)
      detnew[k,i] <- det(m_s2[[i]])
      
      #update distances
      dist[,i] =  sqrt(mahalanobis(x = x, cov = m_s2[[i]], center = v_estmu))
      
      if(detnew[k,i] < det_old){
        det_old <- det(m_s2[[i]])
      }
      else if (detnew[k,i] ==  det_old){
        break
      } 
     }
    lowdet[i] <- detnew[k,i]
  }
  
  #find which initial covariance estimators give us a subset with smallest determinant
  smallest_determinant <- min(lowdet)
  best_start <- which(lowdet == smallest_determinant)
  low <- order(lowdet)[1]

  #raw estimates 
  selected_sample <- x[sort.list(dist[,low])[1:h],]
  raw.center <- apply(selected_sample, 2, mean)
  
  chi_sq <- qchisq (alpha, p)/2
  correction <- alpha / pgamma(chi_sq, p/2 + 1) #fisher consistency correction factor
  raw.cov <- cov(selected_sample)*correction
  
  
  #re weighting
  weights <- rep(0,n)
  dist <- sqrt(mahalanobis(x = x, cov = raw.cov, center = raw.center))
  not_outliers <- vector()
  cut_off <-  sqrt(qchisq(1-delta, df = p))

  for(j in 1:n){
    if(dist[j] < cut_off){
      not_outliers <- c(not_outliers, j)
      weights[j] = 1
    }
  }
  not_outlier_data <- x[weights>0,]
  
  #re weighted estimates 
  center <- apply(not_outlier_data, 2, mean)
  chi_sq <- qchisq (1-delta, p)/2
  correction <- (1-delta) / pgamma(chi_sq, p/2 + 1) #fisher consistency correction factor
  cov <- cov(not_outlier_data)*correction

  #prepare output
  detmcd <- list()
  detmcd$center <- center
  detmcd$cov <- cov
  detmcd$weights <- weights
  detmcd$raw.center <- raw.center
  detmcd$raw.cov <- raw.cov
  detmcd$best <- not_outliers
  detmcd$i_Best <- best_start 
  return(detmcd)
  
  
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

lmDetMCD <- function(x, y, alpha = 0.75, ...) {
  n = length(x)
  data <- cbind(x, y)
  
  #get MCD estimates
  #MCD = covMcd(data, alpha, nsamp = "deterministic")
  MCD = covDetMCD(data, alpha = alpha)
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
