#' Truncated Bernoulli product (TBP) distribution
#'
#' Simulate draws from a TBP distribution
#'
#' @param m number of samples to draw
#' @param pi the TBP parameter vector of length r
#'
#' @returns a matrix of dimension r x m of TBP draws
#' @examples
#' m <- 5000L
#' pi <- c(0.1, 0.05, 0.1)
#' X <- r_tbp(m, pi)
r_tbp <- function(m, pi) {
  Omega <- pmf_tbp(pi)
  s <- sample(x = seq_len(nrow(Omega)), size = m, replace = TRUE, prob = Omega$pmf)
  out <- Omega[s,seq_along(pi)] |> as.matrix() |> t()
  colnames(out) <- NULL
  return(out)
}


#' Returns the pmf of a truncated Bernoulli product (TBP) distribution
#'
#' @param pi parameter of the TBP distribution
#'
#' @returns a data frame with columns (x_1, x_2, \dots, x_r, pmf), where pmf gives the probability of a given (x_1, x_2, \dots, x_r) vector
#' @examples
#' Omega <- pmf_tbp(c(0.1, 0.05, 0.02))
#' Omega <- pmf_tbp(c(0.4, 0.6, 0.2))
pmf_tbp <- function(pi) {
  r <- length(pi)
  Omega <- rep(list(c(0L, 1L)), r) |> setNames(paste0("x_", seq_len(r))) |> expand.grid()
  Omega <- Omega[rowSums(Omega) >= 1L,]
  denom <- 1 - prod(1 - pi)
  pmf <- apply(X = Omega, MARGIN = 1, FUN = function(x) {
    prod(ifelse(x, pi, 1 - pi))
  })/denom
  Omega$pmf <- pmf
  return(Omega)
}


#' Fit truncated Bernoulli product (TBP) model via MLE
#'
#' Fits a TBP model to an r x m matrix of occupancies
#'
#' @param X binary occupancy matrix of dimension r (number of replicates) by m (number of windows)
#'
#' @returns the fitted MLE pi-hat (or a vector containing all NAs if the regularity conditions fail to hold)
#' @examples
#' m <- 5000L
#' pi <- c(0.2, 0.12, 0.04)
#' X <- r_tbp(m, pi)
#' pi_hat <- fit_tpb_model(X)
fit_tpb_model <- function(X) {
  # unconstrained estimate
  q_hat <- rowMeans(X)
  # verify regularity conditions
  if (sum(q_hat) > 1 && all(q_hat > 0) && all(q_hat < 1)) {
    g <- function(c) 1 - prod(1 - c * q_hat) - c
    root_res <- uniroot(g, lower = 1e-10, upper = 1 - 1e-10)
    pi_hat <- root_res$root * q_hat
  } else {
    pi_hat <- rep(NA_real_, nrow(X))
  }
  return(pi_hat)
}


#' shifted negative binomial (SNB) distribution
#'
#' Sample m_plus observations from an SNB distribution with parameters (mu, theta)
#'
#' @param m_plus the number of samples to generate (typically equal to the number of occupied windows)
#' @param mu mean parameter
#' @param theta size parameter
#'
#' @returns a vector of snb variates
#' @examples
#' mu <- 10
#' theta <- 0.5
#' m_plus <- 1000L
#' y_plus <- r_snb(m_plus, mu, theta)
r_snb <- function(m_plus, mu, theta) {
  MASS::rnegbin(n = m_plus, mu = mu, theta = theta) + 1L
}


#' Simulate multi-replicate guide-seq data
#'
#' @param pi parameter vector of TBP model
#' @param mu_vect vector of mu parameters for SNB models
#' @param theta_vect vector of theta parameters for SNB models
#'
#' @returns
#' @examples
#' # NULL DATA
#' pi <- c(0.05, 0.1, 0.02)
#' mu_vect <- c(10, 6, 15)
#' theta_vect <- c(2, 5, 0.3)
#' m <- 10000
#' null_dat <- simulate_multirep_guideseq_data(pi, mu_vect, theta_vect, m)
#'
#' # ALTERNATIVE DATA
#' pi <- c(0.5, 0.6, 0.4)
#' mu_vect <- c(80, 10, 30)
#' theta_vect <- c(20, 21, 15)
#' m <- 15
#' alt_dat <- simulate_multirep_guideseq_data(pi, mu_vect, theta_vect, m)
#'
#' # COMBINED DATA
#' full_dat <- list(X = cbind(alt_dat$X, null_dat$X), Y = cbind( alt_dat$Y, null_dat$Y))
simulate_multirep_guideseq_data <- function(pi, mu_vect, theta_vect, m) {
  X <- r_tbp(m, pi)
  Y <- sapply(X = seq_along(pi), FUN = function(i) {
    y_plus <- r_snb(m_plus = sum(X[i,]), mu = mu_vect[i], theta = theta_vect[i])
    y <- integer(m)
    y[X[i,] == 1L] <- y_plus
    return(y)
  }) |> t()
  out <- list(X = X, Y = Y)
  return(out)
}

