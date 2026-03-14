#' Title
#'
#' @param n_trt
#' @param n_cntrl
#' @param k_trt
#' @param k_cntrl
#' @param rho
#' @param nominal_ci_probability
#'
#' @returns
#' @export
#'
#' @examples
#' p <- 5L
#' amplicon_ids <- factor(x = paste0("amplicon_", seq_len(p)), levels = paste0("amplicon_", seq_len(p)))
#' beta_binom_rho <- rep(5e-4, times = p)
#' data_list <- generate_synthetic_amplicon_seq_data(p = p, r = 3L, pi_cntrl = 0.05, editing_rate = 0.15,
#'                                                   n_amplicons_nonzero_editing = 2L, beta_binom_rho,
#'                                                   amplicon_ids = amplicon_ids)
#' n_trt <- data_list$n_mat_trt[,3]
#' n_cntrl <- data_list$n_mat_cntrl[,3]
#' k_trt <- data_list$k_mat_trt[,3]
#' k_cntrl <- data_list$k_mat_cntrl[,3]
#' rho <- beta_binom_rho[3]
#' alpha_pi <- 1
#' beta_pi <- 100
#' alpha_theta <- 1
#' beta_theta <- 100
compute_bayesian_credible_interval <- function(n_trt, n_cntrl, k_trt, k_cntrl, rho,
                                               alpha_pi, beta_pi, alpha_theta, beta_theta,
                                               nominal_ci_probability = 0.95) {
  # 1. define a grid of points
  u_min <- v_min <- log(10e-5)
  u_max <- v_max <- log(0.999)
  u_grid <- seq(u_min, u_max, length.out = 1000)
  v_grid <- seq(v_min, v_max, length.out = 500)
  exp_u_grid <- exp(u_grid)
  exp_v_grid <- exp(v_grid)

  # 2. compute P1 over the grid
  r <- length(n_trt)
  p1_a <- sapply(X = seq_len(r), FUN = function(i) {
      VGAM::dbetabinom(x = k_trt[i], size = n_trt[i], prob = exp_u_grid, rho = rho, log = TRUE)
    }, simplify = TRUE) |> t() |> colSums()
  p1_b <- alpha_pi * u_grid + (beta_pi - u_grid) * log(1 - exp_u_grid)
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
  theta_posterior <- rowSums(w)
  pi_posterior <- colSums(w)

  # 8. compute the mean, credible interval, and upper credible limit of the posteriors

}
