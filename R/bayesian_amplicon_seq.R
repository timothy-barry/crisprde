#' Compute Bayesian posterior summaries for amplicon-seq editing
#'
#' Approximates the joint posterior of the background mutation rate `pi` and
#' the editing rate `theta` under the beta-binomial model described in the
#' accompanying methods writeup.
#'
#' @param n_trt Numeric vector of total read counts for the treated replicates.
#' @param n_cntrl Numeric vector of total read counts for the control replicates.
#' @param k_trt Numeric vector of mutated read counts for the treated replicates.
#' @param k_cntrl Numeric vector of mutated read counts for the control replicates.
#' @param rho Scalar beta-binomial dispersion parameter shared across replicates.
#' @param alpha_pi First shape parameter of the Beta prior on `pi`.
#' @param beta_pi Second shape parameter of the Beta prior on `pi`.
#' @param alpha_theta First shape parameter of the Beta prior on `theta`.
#' @param beta_theta Second shape parameter of the Beta prior on `theta`.
#' @param alpha Tail probability used for the returned upper credible limits.
#'
#' @returns A list containing posterior means and upper credible limits for
#'   `theta` and `pi`, along with data frames giving the discretized marginal
#'   posterior densities over the evaluation grid.
#' @export
#'
#' @examples
#' data_list <- generate_synthetic_amplicon_seq_data(p = 1L, r = 3L, pi_cntrl = 0.001, editing_rate = 0.025,
#'                                                   n_amplicons_nonzero_editing = 1L, beta_binom_rho = 5e-5,
#'                                                   amplicon_ids = "amplicon_1")
#' n_trt <- data_list$n_mat_trt[,1]
#' n_cntrl <- data_list$n_mat_cntrl[,1]
#' k_trt <- data_list$k_mat_trt[,1]
#' k_cntrl <- data_list$k_mat_cntrl[,1]
#' rho <- 5e-5
#'
#' alpha_pi <- 1
#' beta_pi <- 100
#' alpha_theta <- 1
#' beta_theta <- 1
#'
#' result <- compute_bayesian_credible_interval(n_trt, n_cntrl, k_trt, k_cntrl, rho, alpha_pi, beta_pi, alpha_theta, beta_theta)
#'
compute_bayesian_credible_interval <- function(n_trt, n_cntrl, k_trt, k_cntrl, rho,
                                               alpha_pi, beta_pi, alpha_theta, beta_theta,
                                               alpha = 0.01) {
  # pi -> u -> T
  # theta -> v -> S
  # 1. define a grid of points
  u_min <- v_min <- log(10e-7)
  u_max <- v_max <- log(0.999)
  u_grid <- seq(u_min, u_max, length.out = 500)
  v_grid <- seq(v_min, v_max, length.out = 2000)
  exp_u_grid <- exp(u_grid)
  exp_v_grid <- exp(v_grid)

  # 2. compute P1 over the grid
  r <- length(n_trt)
  p1_a <- sapply(X = seq_len(r), FUN = function(i) {
      VGAM::dbetabinom(x = k_cntrl[i], size = n_cntrl[i], prob = exp_u_grid, rho = rho, log = TRUE)
    }, simplify = TRUE) |> t() |> colSums()
  p1_b <- alpha_pi * u_grid + (beta_pi - 1) * log(1 - exp_u_grid)
  p1 <- p1_a + p1_b

  # 3. compute P2 over grid
  p2 <- alpha_theta * v_grid + (beta_theta - 1) * log(1 - exp_v_grid)

  # 4.compute P1 + P2 outer sum
  p1_p2_outer <- outer(p1, p2, "+")

  # 5. compute outer product over exp(u) + (1 - exp(u)) * exp(v), yielding mean matrix
  mu_mat <- outer(exp_u_grid, exp_v_grid, function(exp_u_grid, exp_v_grid) exp_u_grid + (1 - exp_u_grid) * exp_v_grid)
  p3 <- Reduce(f = "+", x = lapply(X = seq_len(r), FUN = function(i) {
    matrix(data = VGAM::dbetabinom(x = k_trt[i], size = n_trt[i], prob = c(mu_mat), rho = rho, log = TRUE),
           nrow = nrow(mu_mat), ncol = ncol(mu_mat))
  }))

  # 6. compute log-h, max over log h, and the normalized weights
  log_h <- p1_p2_outer + p3
  w <- exp(log_h - max(log_h))
  w <- w/sum(w)

  # 7. compute the approximate posterior density of theta and pi cntrl
  theta_posterior <- colSums(w)
  pi_posterior <- rowSums(w)

  # 8. compute the mean, credible interval, and upper credible limit of the posteriors
  theta_mean <- compute_posterior_mean(theta_posterior, exp_v_grid)
  pi_mean <- compute_posterior_mean(pi_posterior, exp_u_grid)
  theta_upper_credible_limit <- compute_upper_credible_limit(theta_posterior, exp_v_grid, alpha)
  pi_upper_credible_limit <- compute_upper_credible_limit(pi_posterior, exp_u_grid, alpha)
  theta_density_df <- data.frame(theta = exp_v_grid, density = theta_posterior)
  pi_density_df <- data.frame(pi = exp_u_grid, density = pi_posterior)

  # 9. prepare output
  ret <- list(theta_mean = theta_mean,
              theta_upper_credible_limit = theta_upper_credible_limit,
              theta_density_df = theta_density_df,
              pi_mean = pi_mean,
              pi_upper_credible_limit = pi_upper_credible_limit,
              pi_density_df = pi_density_df)
}

compute_posterior_mean <- function(posterior, x_grid) {
  sum(posterior * x_grid)
}

compute_upper_credible_limit <- function(posterior, x_grid, alpha) {
  x_grid[min(which(cumsum(posterior) >= 1 - alpha))]
}
