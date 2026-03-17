#' Run Bayesian amplicon-seq analysis
#'
#' Run Bayesian analysis for each amplicon and return posterior summaries.
#'
#' @param data_list List containing `n_mat_trt`, `n_mat_cntrl`, `k_mat_trt`, `k_mat_cntrl`, and `amplicon_ids`.
#' @param nominal_ci_coverage Credible interval coverage level.
#' @param rho Optional scalar dispersion value to use for every amplicon.
#' @param outlier_mad_thresh MAD threshold used when detecting dispersion outliers.
#' @param min_mutated_read_count Minimum total mutated read count required for an amplicon to contribute to dispersion estimation.
#' @param bias_variance_param Weight on the amplicon-specific pilot dispersion for high-information amplicons.
#' @param alpha_pi,beta_pi Beta prior hyperparameters for `pi`.
#' @param alpha_theta,beta_theta Beta prior hyperparameters for `theta`.
#'
#' @returns A list containing `result_df`, `dispersion_diagnostics`, `theta_posterior_density_df`, and `pi_posterior_density_df`.
#' @export
#'
#' @examples
#' set.seed(2)
#' p <- 25L
#' amplicon_ids <- factor(x = paste0("amplicon_", seq_len(p)), levels = paste0("amplicon_", seq_len(p)))
#' beta_binom_rho <- c(0.005, rep(5e-4, times = p - 1L))
#' data_list <- generate_synthetic_amplicon_seq_data(p = p, r = 3L, pi_cntrl = 0.05, editing_rate = 0.15,
#'                                                   n_amplicons_nonzero_editing = 2L, beta_binom_rho,
#'                                                   amplicon_ids = amplicon_ids)
#'
#' # run analysis
#' alpha_theta <- 1
#' beta_theta <- 1
#' res_bayes <- run_bayesian_amplicon_seq_analysis(data_list = data_list, alpha_theta = alpha_theta, beta_theta = beta_theta)
#' theta_posterior_density_df <- res_bayes$theta_posterior_density_df
#'
#' # make plots
#' make_pilot_dispersion_plot(res_bayes) |> plot()
#' make_amplicon_seq_ci_plot(res_bayes$result_df) |> plot()
#' p_theta_prior <- make_prior_density_plot(alpha = alpha_theta, beta = beta_theta, parameter = "theta")
#' p_theta_posterior_amplicon_1 <- make_posterior_density_plot(posterior_density_df = theta_posterior_density_df, amplicon_id_to_plot = "amplicon_1", parameter = "theta", x_limits = c(0.1, 0.25))
#' p_theta_posterior_amplicon_amplicon_3 <- make_posterior_density_plot(posterior_density_df = theta_posterior_density_df, amplicon_id_to_plot = "amplicon_3", parameter = "theta", x_limits = c(0, 0.02))
run_bayesian_amplicon_seq_analysis <- function(data_list, nominal_ci_coverage = 0.99, rho = NULL,
                                               outlier_mad_thresh = 4, min_mutated_read_count = 50L,
                                               bias_variance_param = 0.5, alpha_pi = 15, beta_pi = 350,
                                               alpha_theta = 1, beta_theta = 1, print_progress = TRUE) {
  # 1. extract the data
  n_mat_trt <- data_list$n_mat_trt
  n_mat_cntrl <- data_list$n_mat_cntrl
  k_mat_trt <- data_list$k_mat_trt
  k_mat_cntrl <- data_list$k_mat_cntrl

  # 2. compute pi tilde, or the Jeffrey-regularized pi estimator (for dispersion estimation)
  estimate_reg_pi_per_amplicon <- function(n_mat, k_mat) (colSums(k_mat) + 0.5)/(colSums(n_mat) + 1)
  pi_tilde_trt_per_amplicon <- estimate_reg_pi_per_amplicon(n_mat = n_mat_trt, k_mat = k_mat_trt)
  pi_tilde_cntrl_per_amplicon <- estimate_reg_pi_per_amplicon(n_mat = n_mat_cntrl, k_mat = k_mat_cntrl)

  # 3. estimate rho via regularized method of moments approach
  dispersion_output <- estimate_dispersions(n_mat_trt = n_mat_trt,
                                            n_mat_cntrl = n_mat_cntrl,
                                            k_mat_trt = k_mat_trt,
                                            k_mat_cntrl = k_mat_cntrl,
                                            pi_tilde_trt_per_amplicon = pi_tilde_trt_per_amplicon,
                                            pi_tilde_cntrl_per_amplicon = pi_tilde_cntrl_per_amplicon,
                                            outlier_mad_thresh = outlier_mad_thresh,
                                            bias_variance_param = bias_variance_param,
                                            rho = rho,
                                            min_mutated_read_count = min_mutated_read_count)
  rho_tilde_per_amplicon <- dispersion_output$rho_tilde_per_amplicon

  # 4. run the Bayesian estimation procedure for each amplicon
  amplicon_wise_res_list <- lapply(X = seq_len(ncol(n_mat_trt)), FUN = function(j) {
    if (print_progress) print(paste0("Running analysis on amplicon ", data_list$amplicon_ids[j]))
    res <- compute_bayesian_estimate_for_amplicon(n_trt = n_mat_trt[,j],
                                                  n_cntrl = n_mat_cntrl[,j],
                                                  k_trt = k_mat_trt[,j],
                                                  k_cntrl = k_mat_cntrl[,j],
                                                  rho = rho_tilde_per_amplicon[j],
                                                  alpha_pi = alpha_pi, beta_pi = beta_pi,
                                                  alpha_theta = alpha_theta, beta_theta = beta_theta,
                                                  alpha = 1 - nominal_ci_coverage)
    amplicon_ids <- data_list$amplicon_ids
    res$amplicon_id <- amplicon_ids[j]
    return(res)
  })

  # 5. prepare the results
  theta_hat <- sapply(X = amplicon_wise_res_list, FUN = function(l) l[["theta_mean"]])
  theta_one_sided_upper_cl <- sapply(X = amplicon_wise_res_list, FUN = function(l) l[["theta_credible_interval"]][["one_sided_bound"]])
  theta_lower_ci <- sapply(X = amplicon_wise_res_list, FUN = function(l) l[["theta_credible_interval"]][["lower_bound"]])
  theta_upper_ci <- sapply(X = amplicon_wise_res_list, FUN = function(l) l[["theta_credible_interval"]][["upper_bound"]])
  theta_posterior_density_df <- lapply(X = amplicon_wise_res_list, FUN = function(l) {
    l[["theta_density_df"]] |> dplyr::mutate(amplicon_id = l[["amplicon_id"]])
  }) |> data.table::rbindlist()
  pi_posterior_density_df <- lapply(X = amplicon_wise_res_list, FUN = function(l) {
    l[["pi_density_df"]] |> dplyr::mutate(amplicon_id = l[["amplicon_id"]])
  }) |> data.table::rbindlist()

  result_df <- data.frame(amplicon_id = data_list$amplicon_ids,
                          rho_hat = rho_tilde_per_amplicon,
                          pilot_rho_hat = dispersion_output$pilot_rho_tilde_per_amplicon,
                          dispersion_outlier = dispersion_output$dispersion_outlier,
                          theta_hat = theta_hat,
                          theta_lower_ci = theta_lower_ci,
                          theta_upper_ci = theta_upper_ci,
                          theta_one_sided_upper_cl = theta_one_sided_upper_cl)
  ret <- list(result_df = result_df,
              dispersion_diagnostics = dispersion_output$dispersion_diagnostics[-1],
              theta_posterior_density_df = theta_posterior_density_df,
              pi_posterior_density_df = pi_posterior_density_df)
  return(ret)
}

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
#' # generate the data
#' data_list <- generate_synthetic_amplicon_seq_data(p = 1L, r = 3L, pi_cntrl = 0.001, editing_rate = 0.0,
#'                                                   n_amplicons_nonzero_editing = 1L, beta_binom_rho = 5e-5,
#'                                                   amplicon_ids = "amplicon_1")
#' n_trt <- data_list$n_mat_trt[,1]
#' n_cntrl <- data_list$n_mat_cntrl[,1]
#' k_trt <- data_list$k_mat_trt[,1]
#' k_cntrl <- data_list$k_mat_cntrl[,1]
#' rho <- 5e-5
#'
#' # set the priors
#' params <- select_beta_hyperparameters(mean = 0.04, variance = 1e-4)
#' alpha_pi <- params[["alpha"]]
#' beta_pi <- params[["beta"]]
#' alpha_theta <- 1
#' beta_theta <- 1
#'
#' # plot the priors
#' p_prior <- cowplot::plot_grid(make_prior_density_plot(alpha = alpha_theta, beta = beta_theta, parameter = "theta"),
#'                               make_prior_density_plot(alpha = alpha_pi, beta = beta_pi, parameter = "pi"))
#' result <- compute_bayesian_estimate_for_amplicon(n_trt, n_cntrl, k_trt, k_cntrl, rho, alpha_pi, beta_pi, alpha_theta, beta_theta)
#' theta_posterior_density_df <- dplyr::mutate(result$theta_density_df, amplicon_id = "amplicon_1")
#' pi_posterior_density_df <- dplyr::mutate(result$pi_density_df, amplicon_id = "amplicon_1")
#' p_posterior <- cowplot::plot_grid(
#'   make_posterior_density_plot(theta_posterior_density_df, amplicon_id_to_plot = "amplicon_1", parameter = "theta", x_limits = c(0, 0.01)),
#'   make_posterior_density_plot(pi_posterior_density_df, amplicon_id_to_plot = "amplicon_1", parameter = "pi", x_limits = c(0, 0.01))
#' )
#'
compute_bayesian_estimate_for_amplicon <- function(n_trt, n_cntrl, k_trt, k_cntrl, rho,
                                                   alpha_pi, beta_pi, alpha_theta, beta_theta,
                                                   alpha = 0.01) {
  # pi -> u -> T
  # theta -> v -> S
  # 1. define a grid of points
  u_min <- v_min <- log(10e-7)
  u_max <- v_max <- log(0.999)
  u_grid <- seq(u_min, u_max, length.out = 500)
  v_grid <- seq(v_min, v_max, length.out = 2500)
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
  theta_credible_interval <- compute_credible_interval(theta_posterior, exp_v_grid, alpha)
  pi_credible_interval <- compute_credible_interval(pi_posterior, exp_u_grid, alpha)
  theta_density_df <- data.frame(theta = exp_v_grid, density = theta_posterior)
  pi_density_df <- data.frame(pi = exp_u_grid, density = pi_posterior)

  # 9. prepare output
  ret <- list(theta_mean = theta_mean,
              theta_credible_interval = theta_credible_interval,
              theta_density_df = theta_density_df,
              pi_mean = pi_mean,
              pi_credible_interval = pi_credible_interval,
              pi_density_df = pi_density_df)
}

compute_posterior_mean <- function(posterior, x_grid) {
  sum(posterior * x_grid)
}

compute_credible_interval <- function(posterior, x_grid, alpha) {
  csum <- cumsum(posterior)
  lower_bound <- x_grid[min(which(csum >= alpha/2))]
  upper_bound <- x_grid[min(which(csum >= 1 - alpha/2))]
  one_sided_bound <- x_grid[min(which(csum >= 1 - alpha))]
  c(lower_bound = lower_bound, upper_bound = upper_bound, one_sided_bound = one_sided_bound)
}

#' Solve for Beta hyperparameters from a target mean and variance
#'
#' Given a desired mean and variance under the Beta parameterization
#' `mean = a / (a + b)` and
#' `variance = ab / ((a + b)^2 (a + b + 1))`,
#' solve for the corresponding positive shape parameters.
#'
#' @param mean Target mean of the Beta distribution. Must lie strictly between 0 and 1.
#' @param variance Target variance of the Beta distribution. Must be positive and strictly smaller than `mean * (1 - mean)`.
#' @returns A named numeric vector with entries `alpha` and `beta`.
#' @export
#'
#' @examples
#' params <- select_beta_hyperparameters(mean = 0.01, variance = 1e-4)
#' p <- make_beta_prior_density_plot(parameter = "pi", alpha = params[["alpha"]], beta = params[["beta"]], xmin = 0, xmax = 0.1)
select_beta_hyperparameters <- function(mean, variance) {
  max_variance <- mean * (1 - mean)
  if (variance >= max_variance) {
    stop("`variance` must be strictly smaller than `mean * (1 - mean)` for a Beta distribution.")
  }
  concentration <- mean * (1 - mean) / variance - 1
  alpha <- mean * concentration
  beta <- (1 - mean) * concentration

  c(alpha = alpha, beta = beta)
}
