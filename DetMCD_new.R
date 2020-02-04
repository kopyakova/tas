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
  index = norms > .Machine$double.eps
  k = z[index,] / norms[index]
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
x <- Eredivisie28
covDetMCD <- function(x, alpha = 0.5) {
  #mahalo distance function
  n = dim(x)[1]
  p = dim(x)[2]
  
  mah_dist_vec <- function(x, m_cov, v_mean){
    v_mah_dist <- vector()
    v_mah_dist <- sqrt(as.matrix(x - v_mean) %*% solve(m_cov) %*% t(as.matrix(x-v_mean)))
    return(v_mah_dist)
  }
  
  qn_estimator <- function(x){
    constant = 2.21914
    n = dim(x)[1]
    p = dim(x)[2]
    h = floor(n/2)+1
    Qn = rep(NA, p)
    for (i in 1:p){
      variable = as.matrix(x[,i])
      rep_matrix = matrix(data = apply(variable, 2, function(x) rep(x, n)), ncol = ncol(variable)*n)
      distance = abs(rep_matrix - t(rep_matrix))
      distance = distance[upper.tri(distance, diag = FALSE)]
      distance = distance[order(distance)]
      
      Qn[i] = constant*distance[h*(h-1)/2];
    }
    return(Qn)
  }
  x <- as.matrix(x)
  v_colmed <- colMedians(x)
  
  qn <- qn_estimator(x)
  m_z <- cbind(
    (x[,1]-v_colmed[1])*(1/qn[1]),
    (x[,2]-v_colmed[2])*(1/qn[2]))
  
  #initial estimates
  m_s <- list()
  m_s[[1]] <- corHT(m_z)
  m_s[[2]] <- corSpearman(m_z)
  m_s[[3]] <- corNSR(m_z)
  m_s[[4]] <- covMSS(m_z) #is wrong
  m_s[[5]] <- covBACON1(m_z)
  m_s[[6]] <- rawCovOGK(m_z)
  # m_z <- list()
  # m_z[[1]] <- m_z_set
  # m_z[[2]] <- m_z_set
  # m_z[[3]] <- m_z_set
  # m_z[[4]] <- m_z_set
  # m_z[[5]] <- m_z_set
  # m_z[[6]] <- m_z_set
  alpha = 0.5
  #get initial values
  
  h <- h.alpha.n(alpha,n,p)
  dist <- matrix(data = 0, ncol = 6, nrow = n)
  d_ik_index <- matrix(data=0, nrow = h, ncol= 6)
  detnew <- matrix(data=0, ncol =6, nrow=50)
  detold <- vector()
  estmu2 <- list()
  m_s2 <- list()
  lowdet <- vector()
  # data = x
  # x = m_z
  for(i in 1:6){
    m_e = eigen(m_s[[i]])$vectors
    m_b = m_z %*% m_e
    m_l = diag(qn_estimator(m_b)^2)# diagonal qn^2 of B
    m_s2[[i]] = m_e %*% m_l %*% t(m_e)
    v_estmu = t(chol2inv(as.matrix(m_s2[[i]])) %*% (colMedians(m_z %*% solve(chol2inv(m_s2[[i]])))))
    for( j in 1:n){
      dist[j,i] = sqrt(mahalanobis(x = m_z[j,], cov = m_s2[[i]], center = v_estmu))
    }
   
    #c-steps

    for(k in 1:25){
      detold[k] <- det(m_s2[[i]])
      selected_sample <- x[order(dist[,i]),]
      selected_sample <- selected_sample[1:h,]
      estmu2[[i]] <- cbind(mean(selected_sample[,1]), mean(selected_sample[,2]))
      v_estmu <- estmu2[[i]] 
      #m_s2[[i]] <- (1/(h-1))*(as.matrix(t(selected_sample-cbind(rep(v_estmu[1], h),rep(v_estmu[2], h)))) %*% as.matrix(selected_sample-cbind(rep(v_estmu[1], h),rep(v_estmu[2], h))))
      m_s2[[i]] <- cov(selected_sample)
      detnew[k,i] <- det(m_s2[[i]])
      for( j in 1:n){
        dist[j,i] = sqrt(mahalanobis(x = x[j,], cov = m_s2[[i]], center = v_estmu))
      }
    
      if(detnew[k,i] == detold[k]){
        break()
      }
    }
    lowdet[i] <- min(detnew[,i][detnew[,i] != 0]) 
  }

  low <- order(lowdet)[1]
  tot = 0
  weight <- vector()
  weights2 <- rep(0,n)
  for(j in 1:n){
    
    if(dist[j,low] <= sqrt(qchisq(.975, df = p))){
      weight <- c(weight, j)
      tot = tot + 1
      weights2[j] = 1
    }
    
  }
  detmcd <- list()
   detmcd$center = cbind(mean(x[weight,1]),mean(x[weight,2]))
   v_estmu = detmcd$center
   selected_sample = x[weight,]
   
   detmcd$cov<- (1/(h-1))*as.matrix(t(selected_sample - rep(v_estmu, each=h))) %*% as.matrix(selected_sample - rep(v_estmu, each=h))
   detmcd$cov<- cov(selected_sample)
   #detmcd$cov = cov(as.matrix(cbind(x[weight,1],x[weight,2])))
   detmcd$weights = weights2
   detmcd$raw.center = estmu2[[low]]
   detmcd$raw.cov = m_s2[[low]]
   detmcd$best = low  
  return(detmcd)
  
  
  # Please note that the subset sizes for the MCD are not simply fractions of 
  # the number of observations in the data set, as discussed in the lectures.
  # You can use function ?h.alpha.n() from package robustbase to compute the 
  # subset size.
  
}
covMcd(x)


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
