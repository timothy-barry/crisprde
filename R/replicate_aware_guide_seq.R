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

# helper function to generate the Omega matrix of combinations
generate_omega <- function(r) {
  Omega <- rep(list(c(0L, 1L)), r) |> setNames(paste0("x_", seq_len(r))) |> expand.grid()
  Omega <- Omega[rowSums(Omega) >= 1L,]
  rownames(Omega) <- NULL
  return(Omega)
}

# helper function to convolve two pmfs
convolve_pmfs <- function(a, b) {
  pmf <- convolve(a, rev(b), type = "open")
  return(pmf)
}

# helper function to convolve a list of pmfs
convolve_pmf_list <- function(pmf_list) {
  if (length(pmf_list) == 1L) {
    out <- pmf_list[[1]]
  } else {
    out <- Reduce(convolve_pmfs, pmf_list)
  }
  out <- pmax(out, 1e-50)
  out <- out / sum(out)
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
  Omega <- generate_omega(r)
  denom <- 1 - prod(1 - pi)
  pmf <- apply(X = Omega, MARGIN = 1, FUN = function(x) {
    prod(ifelse(x, pi, 1 - pi))
  })/denom
  Omega$pmf <- pmf
  return(Omega)
}


#' Fit truncated Bernoulli product (TBP) model via MLE
#'
#' Fits a TBP model to an r x m matrix of occupancies.
#'
#' Returns NULL if the regularity conditions fail to hold.
#'
#' @param X binary occupancy matrix of dimension r (number of replicates) by m (number of windows)
#'
#' @returns the fitted MLE pi-hat (or NULL if the regularity conditions fail to hold)
#' @examples
#' m <- 5000L
#' pi <- c(0.2, 0.12, 0.04)
#' X <- r_tbp(m, pi)
#' pi_hat <- fit_tbp_model(X)
fit_tbp_model <- function(X) {
  # unconstrained estimate
  q_hat <- rowMeans(X)
  # verify regularity conditions
  if (sum(q_hat) > 1 && all(q_hat > 0) && all(q_hat < 1)) {
    g <- function(c) 1 - prod(1 - c * q_hat) - c
    root_res <- uniroot(g, lower = 1e-10, upper = 1 - 1e-10)
    pi_hat <- root_res$root * q_hat
  } else {
    pi_hat <- NULL
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
#' mu_vect <- c(80, 200, 50)
#' theta_vect <- c(20, 21, 15)
#' m <- 15
#' alt_dat <- simulate_multirep_guideseq_data(pi, mu_vect, theta_vect, m)
#'
#' # COMBINED DATA
#' Y_mat <- cbind(alt_dat, null_dat)
#' colnames(Y_mat) <- paste0("window_", seq_len(ncol(Y_mat)))
simulate_multirep_guideseq_data <- function(pi, mu_vect, theta_vect, m) {
  X <- r_tbp(m, pi)
  Y <- sapply(X = seq_along(pi), FUN = function(i) {
    y_plus <- r_snb(m_plus = sum(X[i,]), mu = mu_vect[i], theta = theta_vect[i])
    y <- integer(m)
    y[X[i,] == 1L] <- y_plus
    return(y)
  }) |> t()
  return(Y)
}


#' Run multivariate guide-seq method
#'
#' Runs the multivariate guide-seq method on a replicate-by-window matrix of counts
#'
#' @param Y_mat the r x m integer matrix of UMI counts, where r is the number of replicates and m is the number of windows
#' @param incorporate_occupancy_info a boolean (T/F) indicating whether to incorporate occupancy information into the p-value calculation
#'
#' @returns a data frame containing a p-value for each window
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
#' mu_vect <- c(80, 200, 50)
#' theta_vect <- c(20, 21, 15)
#' m <- 15
#' alt_dat <- simulate_multirep_guideseq_data(pi, mu_vect, theta_vect, m)
#'
#' # COMBINED DATA
#' Y_mat <- cbind(alt_dat, null_dat)
#' incorporate_occupancy_info <- TRUE
run_multivariate_guideseq_method <- function(Y_mat, incorporate_occupancy_info = TRUE, multiplicity_alpha = 0.1) {
  X <- Y_mat > 0
  storage.mode(X) <- "integer"

  # 1. fit occupancy model and compute marginal occupancy probabilities
  if (incorporate_occupancy_info) {
    pi_hat <- fit_tbp_model(X)
    if (!is.null(pi_hat)) tbp_pattern_df <- pmf_tbp(pi_hat)
  }

  # 2. compute shifted NB models and compute convolved pms and cdfs over all combinations
  mu_theta_hat_mat <- apply(X = Y_mat, MARGIN = 1, FUN = function(curr_row) {
    y_plus <- curr_row[curr_row > 0]
    fit <- fit_rob_nb_univariate(y = y_plus - 1)
    fit[c("mu", "theta")]
  }) |> t()
  # generate pmfs
  pmf_list <- apply(X = mu_theta_hat_mat, MARGIN = 1, FUN = function(curr_row) {
    max_count <- qnbinom(p = 1e-50, mu = curr_row[["mu"]], size = curr_row[["theta"]], lower.tail = FALSE)
    dnbinom(x = seq(0L, max_count), mu = curr_row[["mu"]], size = curr_row[["theta"]])
  }, simplify = FALSE)
  # compute all convolution combinations
  Omega <- as.matrix(generate_omega(nrow(Y_mat)))
  conv_pmf_list <- apply(X = Omega, MARGIN = 1, FUN = function(curr_row) {
    convolve_pmf_list(pmf_list[as.logical(curr_row)])
  }, simplify = FALSE)
  # compute right-tail probability list
  right_tail_prob_list <- lapply(X = conv_pmf_list, FUN = function(curr_pmf) {
    rev(cumsum(rev(curr_pmf)))
  })

  # 3. match occupancy pattern to Omega
  omega_keys <- apply(Omega, 1, paste0, collapse = "")
  col_keys <- apply(X, 2, paste0, collapse = "")
  occupancy_pattern_map <- match(col_keys, omega_keys)

  # 3. compute the p-value for each window
  # occupancy + count p-value
  if (incorporate_occupancy_info) {



  } else { # count-only p-value
    p_vals <- sapply(X = seq_len(ncol(Y_mat)), FUN = function(i) {
      y <- Y_mat[,i]
      occupancy_pattern_idx <- occupancy_pattern_map[i]
      n_nonzero <- sum(X[,i])
      test_stat <- sum(y) - n_nonzero
      right_tail_prob_vector <- right_tail_prob_list[[occupancy_pattern_idx]]
      idx <- test_stat + 1L
      if (idx > length(right_tail_prob_vector)) {
        p_val <- right_tail_prob_vector[length(right_tail_prob_vector)]
      } else {
        p_val <- right_tail_prob_vector[idx]
      }
      return(p_val)
    })
  }
  return(p_vals)
}
