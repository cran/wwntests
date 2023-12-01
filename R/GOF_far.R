# V_WS_quantile_far comuputes the 1-alpha qunatile of the beta * chi-squared distribution with nu
#   degrees of freedom, where beta and nu are obtained from a Welch-Satterthwaite approximation
#   of the test statistic V_K. This quantile is used to conduct an approximate size alpha test
#   for the adequacy of FAR(1) models .
# Input: eg = the model residual matrix with functions in columns
#        f = the matrix adjusting the dependence caused by estimating the kernel operator
#        lag = specifies the range of lags 1:K for the test statistic V_K
#        alpha = the significance level to be used in the hypothesis test
#        M = optional argument specifying the sampling size in the related Monte Carlo method
# Output: scalar value of the 1-alpha quantile of the beta * chi-square distribution with nu
#         degrees of freedom (which approximates V_K)
V_WS_quantile_far<-function (eg, f, lag, alpha = 0.05, M = 10000)
{
  mean_V_K <- mean_hat_V_K_far(eg, f, lag)
  var_V_K <- variance_hat_V_K_far(eg, f, lag, M = M)
  beta <- var_V_K/(2 * mean_V_K)
  nu <- 2 * (mean_V_K^2)/var_V_K
  quantile <- beta * qchisq(1 - alpha, nu)
  statistic <- t_statistic_V(eg, lag)
  p_val <- pchisq(statistic/beta, nu, lower.tail = FALSE)
  list(statistic = statistic, quantile = quantile, p_value = p_val)
}

# mean_hat_V_K_far computes the approximation of the mean which is used in the Welch-
#   Satterthwaite approximation as mean of the chi-squared random variable approximating V_K.
# Input: eg = the model residual matrix with functions in columns
#        f = the matrix adjusting the dependence caused by estimating the kernel operator
#        lag = specifies the range of lags 1:K for for the test statistic V_K
# Output: scalar approximation of the mean of the test statistic V_K.
mean_hat_V_K_far<-function (eg, f, lag)
{
  J <- NROW(eg)
  sum1 <- 0
  store <- covariance_diag_store_far(eg, f, lag)
  for (i in 1:lag) {
    sum1 <- sum1 + sum(store[[i]])
  }
  mu_hat_V_K <- (1/(J^2)) * sum1
  mu_hat_V_K
}

# covariance_diag_store_far returns a list storage of the approximate covariances c^hat_i_i(t,s,t,s),
#   for all i in 1:K, for each encoding all values of t,s in U_J X U_J.
# Input: eg = the model residual matrix with functions in columns
#        f = the matrix adjusting the dependence caused by estimating the kernel operator
#        lag = specifies the range of lags 1:K for for the test statistic V_K
# Output: a list containing K 2-D arrays encoding c^hat_i_j(t,s,t,s) evaluated at all (t,s) in
#         U_JxU_J, for i in 1:K
covariance_diag_store_far<-function (eg, f, lag)
{
  cov_i_store <- list()
  for (j in 1:lag) {
    cov_i_store[[j]] <- diagonal_covariance_i_far(eg, f, j)
  }
  cov_i_store
}


# diagonal_covariance_i_j returns the approximate covariance c^hat_i_i(t,s,t,s), encoding all
#   values of t,s in U_J X U_J, i in 1:T.
# Input: eg = the model residual matrix with functions in columns
#        f = the matrix adjusting the dependence caused by estimating the kernel operator
#        lag = specifies the range of lags 1:K for for the test statistic V_K
# Output: a 2-D array encoding c^hat_i_j(t,s,t,s) evaluated at all (t,s) in U_JxU_J.

diagonal_covariance_i_far<-function (eg, f, lag)
{
  N = NCOL(eg)
  J = NROW(eg)

  e_sq_times_e_sq<-(eg[,(1+lag):N])^2%*%t((eg[,1:(N-lag)])^2)/N
  e_sq_times_f_sq<-(eg[,(1+lag):N])^2%*%(f[lag:(N-1),,lag])^2/N
  e_sq_ef<- (eg[,(1+lag):N])^2%*%t(eg[,1:(N-lag)]*t(f[lag:(N-1),,lag]))/N

  ee_ef_sq<-eg[,(1+lag):N]%*%t(eg[,1:(N-lag)])/N - eg[,(1+lag):N]%*%f[lag:(N-1),,lag]/N

  cov<-e_sq_times_e_sq + e_sq_times_f_sq - 2*e_sq_ef - (ee_ef_sq)^2
  cov
}


# variance_hat_V_K_far computes the approximation of the variance which is used in
#   the Welch- Satterthwaite approximation as the variance of the chi-squared random variable
#   approximating V_K.
# Input: eg = the model residual matrix with functions in columns
#        f = the matrix adjusting the dependence caused by estimating the kernel operator
#        lag = specifies the range of lags 1:K for for the test statistic V_K
#        M = optional argument specifying the sampling size in the related Monte Carlo method
# Output: scalar approximation of the variance of the test statistic V_K.
variance_hat_V_K_far<-function (eg, f, lag, M = NULL)
{
  N <- NCOL(eg)
  K=lag
  sum1 <- 0
  for (i in 1:K) {
    sum1 <- sum1 + MCint_eta_approx_i_j_far(eg, f, i, i, M = M)
  }
  bandwidth <- ceiling(0.25 * (N^(1/3)))
  if (K > 1) {
    for (i in 1:(K - 1)) {
      for (j in (i + 1):K) {
        if (abs(i - j) > bandwidth) {
          next
        }
        sum1 <- sum1 + (2 * MCint_eta_approx_i_j_far(eg, f, i, j, M = M))
      }
    }
  }
  variance_V_K <- sum1
  variance_V_K
}


# MCint_eta_approx_i_j_far computes an approximation  using the
#   Monte Carlo integration method "MCint"..
# Input: eg = the model residual matrix with functions in columns
#        f = the matrix adjusting the dependence caused by estimating the kernel operator
#        lag = specifies the range of lags 1:K for for the test statistic V_K
#        i,j = the indices i,j in 1:T that we are computing eta^hat_i_j for
#        M = number of vectors (v1, v2, v3, v4) to sample uniformly from U_J X U_J X U_J X U_J
# Output: scalar value of eta^_hat_i_j computed using the MCint method.
MCint_eta_approx_i_j_far<-function (eg, f, i, j, M = NULL)
{
  J <- NROW(eg)
  N <- NCOL(eg)
  if (is.null(M)) {
    M = floor((max(150 - N, 0) + max(100 - J, 0) + (J/sqrt(2))))
  }

  rand_samp_mat <- matrix(nrow = M, ncol = 4)
  rand_samp_mat <- cbind(sample(1:J, M, replace = TRUE),sample(1:J, M, replace = TRUE),sample(1:J, M, replace = TRUE),sample(1:J, M, replace = TRUE))

  eta_hat_i_j_sum <- 0
  for (k in 1:M) {
    cov <- scalar_covariance_i_j_far(eg, f, i, j, rand_samp_mat[k, ])
    eta_hat_i_j_sum <- eta_hat_i_j_sum + (cov^2)
  }
  eta_hat_i_j <- (2/M) * eta_hat_i_j_sum
  eta_hat_i_j
}



# scalar_covariance_i_j_far returns the approximate covariance c^hat_i_j(t,s,u,v) evaluated at a
#   given t,s,u,v in U_J X U_J X U_J X U_J (for use in MCint method).
# Input: eg = the model residual matrix with functions in columns
#        f = the matrix adjusting the dependence caused by estimating the kernel operator
#        i,j = the indices i,j in 1:N that we are computing the covariance for
#        times = a 4-element vector representing the values (t,s,u,v)
# Output: scalar value of the computed covariance c^hat_i_j(t,s,u,v).
#
scalar_covariance_i_j_far<-function (eg, f, i, j, times)
{
  J <- NROW(eg)
  N <- NCOL(eg)

  k=1+max(i,j)

  iuv<-eg[times[1], k:N] * eg[times[2], (k-i):(N-i)] - eg[times[1], k:N] * f[(k-1):(N-1), times[2], i]

  juv<-eg[times[3], k:N] * eg[times[4], (k-j):(N-j)] - eg[times[3], k:N] * f[(k-1):(N-1), times[4], j]

  Eiuv<-t(eg[times[1],(1+i):N])%*%eg[times[2],1:(N-i)]/N - t(eg[times[1], (1+i):N])%*%f[i:(N-1), times[2], i]/N

  Ejuv<-t(eg[times[3],(1+j):N])%*%eg[times[4],1:(N-j)]/N - t(eg[times[3], (1+j):N])%*%f[j:(N-1), times[4], j]/N

  cov <- sum(iuv * juv)/N - Eiuv*Ejuv
  cov
}

