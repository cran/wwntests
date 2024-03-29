# V_WS_hyp_test comuputes the 1-alpha qunatile of the beta * chi-squared distribution with nu
#   degrees of freedom, where beta and nu are obtained from a Welch-Satterthwaite approximation
#   of the test statistic V_K. This quantile is used to conduct an approximate size alpha test
#   of the hypothesis H'_0_K.
# Input: f_data = the functional data matrix with functions in columns
#        K = specifies the range of lags 1:K for the test statistic V_K
#        alpha = the significance level to be used in the hypothesis test
#        M = optional argument specifying the sampling size in the related Monte Carlo method
# Output: scalar value of the 1-alpha quantile of the beta * chi-square distribution with nu
#         degrees of freedom (which approximates V_K)
V_WS_quantile <- function(f_data, K, alpha=0.05, M=NULL) {
  mean_V_K <- mean_hat_V_K(f_data, K)
  var_V_K <- variance_hat_V_K(f_data, K, M=M)
  beta <- var_V_K / (2 * mean_V_K)
  nu <- 2 * (mean_V_K^2) / var_V_K
  quantile <- beta * qchisq(1 - alpha, nu)
  statistic <- t_statistic_V(f_data, K)
  p_val <- pchisq(statistic / beta, nu, lower.tail = FALSE)
  list(statistic = statistic, quantile = quantile, p_value = p_val)
}

V_WS_quantile_iid <- function(f_data, K, alpha=0.05) {
  mean_V_K <- mean_hat_V_K_iid(f_data, K)
  var_V_K <- variance_hat_V_K_iid(f_data, K)
  beta <- var_V_K / (2 * mean_V_K)
  nu <- 2 * (mean_V_K^2) / var_V_K
  quantile <- beta * qchisq(1 - alpha, nu)
  statistic <- t_statistic_V(f_data, K)
  p_val <- pchisq(statistic / beta, nu, lower.tail = FALSE)
  list(statistic = statistic, quantile = quantile, p_value = p_val)
}

# Q_WS_hyp_test comuputes the 1-alpha qunatile of the beta * chi-squared distribution with nu
#   degrees of freedom, where beta and nu are obtained from a Welch-Satterthwaite approximation
#   of the test statistic Q_h. This quantile is used to conduct an approximate size alpha test
#   of the hypothesis H_0_h.
# Input: f_data = the functional data matrix with functions in columns
#        lag = specifies the lag used for the test statistic Q_h
#        alpha = the significance level to be used in the hypothesis test
#        M = optional argument specifying the sampling size in the related Monte Carlo method
# Output: scalar value of the 1-alpha quantile of the beta * chi-square distribution with nu
#         degrees of freedom (which approximates Q_h).
Q_WS_quantile <- function(f_data, lag, alpha=0.05, M=NULL) {
  mean_Q_h <- mean_hat_Q_h(f_data, lag)
  var_Q_h <- variance_hat_Q_h(f_data, lag, M=M)
  beta <- var_Q_h / (2 * mean_Q_h)
  nu <- 2 * (mean_Q_h^2) / var_Q_h
  quantile <- beta * qchisq(1 - alpha, nu)
  statistic <- t_statistic_Q(f_data, lag)
  p_val <- pchisq(statistic / beta, nu, lower.tail = FALSE)
  list(statistic = statistic, quantile = quantile, p_value = p_val)
}

# Q_WS_quantile_iid computes the size alpha test of the hypothesis H_0_h using the WS
#   Approximation under the assumption that the data follows a strong white noise.
# Input: f_data = the functional data matrix with functions in columns
#        alpha = the significance level to be used in the hypothesis test
# Output: scalar value of the 1-alpha quantile of the beta * chi-square distribution with nu
#         degrees of freedom (which approximates Q_h) (computed under a strong white noise
#         assumption).
Q_WS_quantile_iid <- function(f_data, alpha=0.05) {
  mean_Q_h <- mean_hat_Q_h_iid(f_data)
  var_Q_h <- variance_hat_Q_h_iid(f_data)
  beta <- var_Q_h / (2 * mean_Q_h)
  nu <- 2 * (mean_Q_h^2) / var_Q_h
  quantile <- beta * qchisq(1 - alpha, nu)
  statistic <- t_statistic_Q(f_data, lag = 1)
  p_val <- pchisq(statistic / beta, nu, lower.tail = FALSE)
  list(statistic = statistic, quantile = quantile, p_value = p_val)
}


#' Compute size alpha single-lag hypothesis test under weak or strong white noise assumption
#'
#' `Q_WS_hyp_test` computes the size alpha test of a single lag hypothesis under a weak white noise
#' or strong white noise assumption using a Welch-Satterthwaite Approximation.
#'
#' @param f_data The functional data matrix with observed functions in the columns
#' @param lag Positive integer value. The lag to use to compute the single lag test statistic.
#' @param alpha Numeric value between 0 and 1 specifying the significance level to be used in the specified
#' hypothesis test. The default value is 0.05. Note, the significance value is only ever used to compute the
#' 1-alpha quantile of the limiting distribution of the specified test's test statistic.
#' @param iid A Boolean value, FALSE by default. If given TRUE, the hypothesis test will use a strong-white
#' noise assumption (instead of a weak-white noise assumption).
#' @param M Positive integer value. Number of Monte-Carlo simulations for the Welch-Satterthwaite approximation.
#' @param bootstrap A Boolean value, FALSE by default. If given TRUE, the hypothesis test is done by
#' approximating the limiting distribution of the test statistic via a block bootstrap process.
#' @param block_size A positive Integer value, with the default value being computed via the adaptive
#' bandwidth selection method in the "spectral" test. Determines the block size (of each block in each
#' bootstrap sample) if the test is being bootstrapped.
#' @param straps A positive Integer, with a default value of 300. Determines the number of bootstrap samples
#' to take if the test is being bootstrapped. Only used if 'bootstrap' == TRUE.
#' @param moving A Boolean value, FALSE by default. If given TRUE, the performed block bootstrap will be moving
#' rather than stationary.
#' @return A list containing the p-value, the quantile, and a boolean value indicating whether or not the
#' hypothesis is rejected.
#'
#' @import stats
Q_WS_hyp_test <- function(f_data, lag, alpha=0.05, iid=FALSE,
                          M=NULL, bootstrap=FALSE,
                          block_size='adaptive', straps=300, moving = FALSE) {
  statistic <- t_statistic_Q(f_data, lag)
  if (bootstrap == TRUE) {
    if (block_size == 'adaptive') {
      block_size <- ceiling(adaptive_bandwidth(f_data, kernel = 'Bartlett'))
    }
    bootsraps <- list()
    bootstrap_samples <- block_bootsrap(f_data, block_size, B = straps, moving = moving)
    stats_distr <- lapply(bootstrap_samples, t_statistic_Q, lag=lag)
    statistic <- t_statistic_Q(f_data, lag=lag)
    quantile <- quantile(as.numeric(stats_distr), 1 - alpha)
    p_value <- sum(statistic > stats_distr) / length(stats_distr)
    list(statistic = as.numeric(statistic), quantile = as.numeric(quantile),
         p_value = as.numeric(p_value), block_size = block_size)
  } else if (iid == FALSE) {
    results <- Q_WS_quantile(f_data, lag, alpha=alpha, M=M)
    statistic <- results$statistic
    quantile <- results$quantile
    p_val <- results$p_val
    reject <- statistic > quantile
    list(statistic = statistic, quantile = quantile, p_value = p_val)
  } else {
    results <- Q_WS_quantile_iid(f_data, alpha=alpha)
    statistic <- results$statistic
    quantile <- results$quantile
    p_val <- results$p_val
    reject <- statistic > quantile
    list(statistic= statistic, quantile = quantile, p_value = p_val)
  }
}
