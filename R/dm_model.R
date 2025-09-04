#' Fit DM model
#'
#' Fits a DM model to a matrix of counts X.
#'
#' @param X A matrix of dimension n by k, where n is the number of samples and k the number of categories
#'
#' @returns a fitted DM model
#' @export
#'
#' @examples
#' n_umis_in_window <- 500
#' n <- 1000
#' K <- 5
#' pi <- c(0.1, 0.2, 0.4, 0.2, 0.1)
#' X <- dirmult::simPop(J = n, K = K, n = n_umis_in_window, pi = pi, theta = 0.05)$data
#' fit <- fit_dm_model(X)
#' pi_hat <- fit$pi
#' theta_hat <- fit$theta
#' p_vals <- compute_lrt_p_vals(X = X, pi = pi_hat, theta = theta_hat)
fit_dm_model <- function(X) {
  fit <- dirmult::dirmult(data = X, trace = FALSE)
  return(fit)
}

#' Compute DM log likelihood
#'
#' Computes the DM log-likelihood for a single vector of counts
#'
#' @param r the vector of counts
#' @param pi the mean of the DM model
#' @param theta the size parameter of the DM model
#'
#' @returns the log-likelihood (up to the multinomial coefficient normalizing constant)
compute_dm_log_lik <- function(r, pi, theta) {
  total <- sum(r)
  gamma_plus <- (1 - theta)/theta
  gamma <- gamma_plus * pi
  log_lik <- lgamma(gamma_plus) - lgamma(total + gamma_plus) + sum(lgamma(r + gamma) - lgamma(gamma))
  return(log_lik)
}

compute_lrt_p_vals <- function(X, pi, theta) {
  K <- ncol(X)
  phi <- (1 - theta)/theta
  p_values <- apply(X, MARGIN = 1, FUN = function(r) {
    reduced_model_log_lik <- compute_dm_log_lik(r, pi, theta)
    pi_hat <- (r + phi/K)/(sum(r) + phi)
    full_model_log_lik <- compute_dm_log_lik(r, pi_hat, theta)
    chi_sq <- -2 * (reduced_model_log_lik - full_model_log_lik)
    pchisq(q = chi_sq, df = K - 1, lower.tail = FALSE)
  })
  return(p_values)
}
