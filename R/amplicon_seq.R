#' Run amplicon-seq analysis
#'
#' Compute Wald-style inference for each amplicon under the beta-binomial model.
#' Dispersion is estimated only from amplicons with at least
#' `min_mutated_read_count` mutated reads across treatment and control
#' replicates; lower-information amplicons are retained in the output and are
#' assigned the shared dispersion estimate.
#'
#' @param data_list A list containing `n_mat_trt`, `n_mat_cntrl`, `k_mat_trt`, `k_mat_cntrl`, and `amplicon_ids`.
#' @param editing_threshold Scalar in `[0, 1]`; test whether the editing rate is greater than this threshold.
#' @param nominal_ci_coverage Nominal coverage level for the reported confidence intervals.
#' @param nominal_fdr Nominal false discovery rate for the BH-adjusted significance calls.
#' @param rho Optional scalar dispersion parameter. If supplied, this value is used for every amplicon instead of estimating dispersion from the data.
#' @param tail Either `"right"` or `"both"`, indicating the tail used to compute p-values.
#' @param outlier_mad_thresh MAD-based threshold used to flag unusually large pilot dispersion estimates before computing the shared dispersion estimate.
#' @param min_mutated_read_count Minimum total mutated read count required for an amplicon to contribute to dispersion estimation.
#' @param bias_variance_param Scalar in `[0, 1]` giving the weight on the pilot dispersion estimate for high-information amplicons. A value of 0 uses only the shared dispersion estimate, while 1 uses only the pilot estimate.
#' @returns A list with components `result_df` and `dispersion_diagnostics`.
#'   `result_df` contains one row per amplicon with the estimated dispersion,
#'   clipped editing-rate estimate, confidence interval, p-value, and
#'   significance call. `dispersion_diagnostics` summarizes the pilot
#'   dispersion distribution and the shared dispersion estimate when `rho` is
#'   estimated from the data.
#' @export
#'
#' @examples
#' # generate data
#' set.seed(2)
#' p <- 25L
#' amplicon_ids <- factor(x = paste0("amplicon_", seq_len(p)), levels = paste0("amplicon_", seq_len(p)))
#' beta_binom_rho <- c(0.005, rep(5e-4, times = p - 1L))
#' data_list <- generate_synthetic_amplicon_seq_data(p = p, r = 3L, pi_cntrl = 0.05, editing_rate = 0.15,
#'                                                   n_amplicons_nonzero_editing = 2L, beta_binom_rho,
#'                                                   amplicon_ids = amplicon_ids)
#'
#' # run analysis
#' res_freq <- run_freqentist_amplicon_seq_analysis(data_list)
#'
#' # make plots
#' make_pilot_dispersion_plot(res_freq) |> plot()
#' make_amplicon_seq_ci_plot(res_freq$result_df) |> plot()
#' make_amplicon_seq_p_value_plot(res_freq$result_df) |> plot()
run_freqentist_amplicon_seq_analysis <- function(data_list, editing_threshold = 0, nominal_ci_coverage = 0.99,
                                                 nominal_fdr = 0.1, rho = NULL, tail = "right",
                                                 outlier_mad_thresh = 4, min_mutated_read_count = 50L,
                                                 bias_variance_param = 0.5) {
  # 1. extract the data
  n_mat_trt <- data_list$n_mat_trt
  n_mat_cntrl <- data_list$n_mat_cntrl
  k_mat_trt <- data_list$k_mat_trt
  k_mat_cntrl <- data_list$k_mat_cntrl

  # 2. estimate of pi and theta
  estimate_pi_per_amplicon <- function(n_mat, k_mat) colSums(k_mat)/colSums(n_mat)
  pi_hat_trt_per_amplicon <- estimate_pi_per_amplicon(n_mat = n_mat_trt, k_mat = k_mat_trt)
  pi_hat_cntrl_per_amplicon <- estimate_pi_per_amplicon(n_mat = n_mat_cntrl, k_mat = k_mat_cntrl)
  theta_hat_per_amplicon <- (pi_hat_trt_per_amplicon - pi_hat_cntrl_per_amplicon)/(1 - pi_hat_cntrl_per_amplicon)

  # 3. compute pi tilde, or the Jeffrey-regularized pi estimator
  estimate_reg_pi_per_amplicon <- function(n_mat, k_mat) (colSums(k_mat) + 0.5)/(colSums(n_mat) + 1)
  pi_tilde_trt_per_amplicon <- estimate_reg_pi_per_amplicon(n_mat = n_mat_trt, k_mat = k_mat_trt)
  pi_tilde_cntrl_per_amplicon <- estimate_reg_pi_per_amplicon(n_mat = n_mat_cntrl, k_mat = k_mat_cntrl)

  # 4. estimate rho via method of moments estimator
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
  pilot_rho_tilde_per_amplicon <- dispersion_output$pilot_rho_tilde_per_amplicon
  dispersion_outlier <- dispersion_output$dispersion_outlier
  dispersion_diagnostics <- dispersion_output$dispersion_diagnostics
  rho_tilde_per_amplicon <- dispersion_output$rho_tilde_per_amplicon

  # 5. compute the standard error of the pi_hats
  compute_pi_tilde_ses <- function(n_mat, pi_tilde_per_amplicon, rho_tilde_per_amplicon) {
    sapply(X = seq_len(ncol(n_mat)), FUN = function(i) {
      n <- n_mat[,i]
      n_tot <- sum(n)
      pi_tilde <- pi_tilde_per_amplicon[[i]]
      rho_tilde <- rho_tilde_per_amplicon[i]
      sqrt(sum(n * pi_tilde * (1 - pi_tilde) * (1 + (n - 1) * rho_tilde)))/n_tot
    })
  }
  pi_tilde_se_trt_per_amplicon <- compute_pi_tilde_ses(n_mat = n_mat_trt,
                                                       pi_tilde_per_amplicon = pi_tilde_trt_per_amplicon,
                                                       rho_tilde_per_amplicon = rho_tilde_per_amplicon)
  pi_tilde_se_cntrl_per_amplicon <- compute_pi_tilde_ses(n_mat = n_mat_cntrl,
                                                         pi_tilde_per_amplicon = pi_tilde_cntrl_per_amplicon,
                                                         rho_tilde_per_amplicon = rho_tilde_per_amplicon)

  # 6. compute the theta_hat standard errors
  theta_tilde_se_per_amplicon <- sapply(X = seq_along(pi_tilde_trt_per_amplicon), FUN = function(i) {
    pi_tilde_trt <- pi_tilde_trt_per_amplicon[i]
    pi_tilde_cntrl <- pi_tilde_cntrl_per_amplicon[i]
    pi_tilde_se_trt <- pi_tilde_se_trt_per_amplicon[i]
    pi_tilde_se_cntrl <- pi_tilde_se_cntrl_per_amplicon[i]
    1/(1 - pi_tilde_cntrl)^2 * sqrt((1 - pi_tilde_cntrl)^2 * pi_tilde_se_trt^2 + (1 - pi_tilde_trt)^2 * pi_tilde_se_cntrl^2)
  })

  # 7. compute confidence intervals and standard errors
  mult_factor_two_sided <- qnorm(p = (1 - nominal_ci_coverage)/2, lower.tail = FALSE)
  mult_factor_one_sided <- qnorm(p = 1 - nominal_ci_coverage, lower.tail = FALSE)
  lower_theta_ci_per_amplicon <- pmax(pmin(theta_hat_per_amplicon - mult_factor_two_sided * theta_tilde_se_per_amplicon, 1), 0)
  upper_theta_ci_per_amplicon <- pmax(pmin(theta_hat_per_amplicon + mult_factor_two_sided * theta_tilde_se_per_amplicon, 1), 0)
  one_sided_upper_cl_per_amplicon <- pmax(pmin(theta_hat_per_amplicon + mult_factor_one_sided * theta_tilde_se_per_amplicon, 1), 0)
  z_score_per_amplicon <- (theta_hat_per_amplicon - editing_threshold)/theta_tilde_se_per_amplicon
  p_val_per_amplicon <- if (tail == "right") {
    pnorm(q = z_score_per_amplicon, lower.tail = FALSE)
  } else if (tail == "both") {
    2 * pmin(pnorm(q = z_score_per_amplicon, lower.tail = FALSE), pnorm(q = z_score_per_amplicon, lower.tail = TRUE))
  } else {
    stop("Tail not recognized.")
  }
  significant <- p.adjust(p = p_val_per_amplicon, method = "BH") <= nominal_fdr

  # 8. return output
  to_return <- data.frame(amplicon_id = data_list$amplicon_ids,
                          rho_hat = rho_tilde_per_amplicon,
                          pilot_rho_hat = pilot_rho_tilde_per_amplicon,
                          dispersion_outlier = dispersion_outlier,
                          theta_hat = pmax(pmin(theta_hat_per_amplicon, 1), 0),
                          theta_lower_ci = lower_theta_ci_per_amplicon,
                          theta_upper_ci = upper_theta_ci_per_amplicon,
                          theta_one_sided_upper_cl = one_sided_upper_cl_per_amplicon,
                          pi_hat_cnrl = pi_hat_cntrl_per_amplicon,
                          p_value = pmax(p_val_per_amplicon, 1e-250),
                          significant = significant)
    out <- list(result_df = to_return,
                dispersion_diagnostics = dispersion_diagnostics[-1])
    return(out)
}


estimate_dispersions <- function(n_mat_trt, n_mat_cntrl, k_mat_trt, k_mat_cntrl,
                                 pi_tilde_trt_per_amplicon, pi_tilde_cntrl_per_amplicon,
                                 outlier_mad_thresh, bias_variance_param,
                                 rho, min_mutated_read_count) {
  pilot_rho_tilde_per_amplicon <- rep(NA_real_, length = ncol(n_mat_trt))
  dispersion_outlier <- rep(NA, length = ncol(n_mat_trt))
  dispersion_diagnostics <- NULL

  #  find high-information amplicons
  high_info_content <- which(colSums(k_mat_trt) + colSums(k_mat_cntrl) >= min_mutated_read_count) |> as.integer()
  if (is.null(rho) && length(high_info_content) == 0L) {
    stop("No amplicon has sufficiently many mutated reads to estimate a dispersion parameter. Consider specifying rho via the `rho` argument.")
  }

  if (is.null(rho)) {
    # get the pilot estimate of rho for each amplicon with high information content
    pilot_rho_tilde_per_amplicon_high_info <- estimate_rho_per_amplicon(n_mat_trt = n_mat_trt[,high_info_content, drop = FALSE],
                                                                        n_mat_cntrl = n_mat_cntrl[,high_info_content, drop = FALSE],
                                                                        k_mat_trt = k_mat_trt[,high_info_content, drop = FALSE],
                                                                        k_mat_cntrl = k_mat_cntrl[,high_info_content, drop = FALSE],
                                                                        pi_tilde_trt_per_amplicon = pi_tilde_trt_per_amplicon[high_info_content],
                                                                        pi_tilde_cntrl_per_amplicon = pi_tilde_cntrl_per_amplicon[high_info_content])
    dispersion_diagnostics <- flag_outlier_dispersions(pilot_rho_tilde_per_amplicon = pilot_rho_tilde_per_amplicon_high_info,
                                                       outlier_mad_thresh = outlier_mad_thresh)
    disp_ok_v <- dispersion_diagnostics$ok_disp
    shared_rho_tilde <- estimate_global_rho(n_mat_trt = n_mat_trt[,high_info_content[disp_ok_v], drop = FALSE],
                                            n_mat_cntrl = n_mat_cntrl[,high_info_content[disp_ok_v], drop = FALSE],
                                            k_mat_trt = k_mat_trt[,high_info_content[disp_ok_v], drop = FALSE],
                                            k_mat_cntrl = k_mat_cntrl[,high_info_content[disp_ok_v], drop = FALSE],
                                            pi_tilde_trt_per_amplicon = pi_tilde_trt_per_amplicon[high_info_content[disp_ok_v], drop = FALSE],
                                            pi_tilde_cntrl_per_amplicon = pi_tilde_cntrl_per_amplicon[high_info_content[disp_ok_v], drop = FALSE])
    dispersion_diagnostics$shared_rho_hat <- shared_rho_tilde
    weight_vector <- ifelse(disp_ok_v, bias_variance_param, 1)
    rho_tilde_per_amplicon_high_info <- weight_vector * pilot_rho_tilde_per_amplicon_high_info + (1 - weight_vector) * shared_rho_tilde

    # expand results to entire set of amplicons
    rho_tilde_per_amplicon <- rep(shared_rho_tilde, length = ncol(n_mat_trt))
    rho_tilde_per_amplicon[high_info_content] <- rho_tilde_per_amplicon_high_info
    dispersion_outlier[high_info_content] <- !disp_ok_v
    pilot_rho_tilde_per_amplicon[high_info_content] <- pilot_rho_tilde_per_amplicon_high_info
  } else {
    rho_tilde_per_amplicon <- rep(rho, ncol(n_mat_trt))
  }
  return(list(pilot_rho_tilde_per_amplicon = pilot_rho_tilde_per_amplicon,
              dispersion_outlier = dispersion_outlier,
              dispersion_diagnostics = dispersion_diagnostics,
              rho_tilde_per_amplicon = rho_tilde_per_amplicon))
}


estimate_global_rho <- function(n_mat_trt, n_mat_cntrl, k_mat_trt, k_mat_cntrl,
                                pi_tilde_trt_per_amplicon, pi_tilde_cntrl_per_amplicon) {
  n_vect <- c(as.integer(n_mat_trt), as.integer(n_mat_cntrl))
  k_vect <- c(as.integer(k_mat_trt), as.integer(k_mat_cntrl))
  p <- ncol(n_mat_trt)
  r <- nrow(n_mat_trt)
  pi_tilde_vect <- c(rep(pi_tilde_trt_per_amplicon, each = r), rep(pi_tilde_cntrl_per_amplicon, each = r))
  gamma <- 2 * p * (r - 1)
  shifted_pearson_stat <- function(cand_rho, n_vect, k_vect, pi_tilde_vect, gamma) {
    sum((k_vect - n_vect * pi_tilde_vect)^2/(n_vect * pi_tilde_vect * (1 - pi_tilde_vect) * (1 + (n_vect - 1) *  cand_rho))) - gamma
  }
  rho_tilde <- if (shifted_pearson_stat(0, n_vect, k_vect, pi_tilde_vect, gamma) <= 0) {
    0
  } else {
    uniroot(f = shifted_pearson_stat, interval = c(0, 0.5), n_vect, k_vect, pi_tilde_vect, gamma)$root
  }
  return(rho_tilde)
}


estimate_rho_per_amplicon <- function(n_mat_trt, n_mat_cntrl, k_mat_trt, k_mat_cntrl,
                                      pi_tilde_trt_per_amplicon, pi_tilde_cntrl_per_amplicon) {
  p <- ncol(n_mat_trt)
  r <- nrow(n_mat_trt)
  sapply(X = seq_len(p), FUN = function(amplicon_idx) {
    n_vect <- c(n_mat_trt[,amplicon_idx], n_mat_cntrl[,amplicon_idx])
    k_vect <- c(k_mat_trt[,amplicon_idx], k_mat_cntrl[,amplicon_idx])
    gamma <- 2 * r - 2
    pi_tilde_vect <- c(rep(pi_tilde_trt_per_amplicon[amplicon_idx], each = r),
                       rep(pi_tilde_cntrl_per_amplicon[amplicon_idx], each = r))
    shifted_pearson_stat <- function(cand_rho, n_vect, k_vect, pi_tilde_vect, gamma) {
      sum((k_vect - n_vect * pi_tilde_vect)^2/(n_vect * pi_tilde_vect * (1 - pi_tilde_vect) * (1 + (n_vect - 1) *  cand_rho))) - gamma
    }
    if (shifted_pearson_stat(0, n_vect, k_vect, pi_tilde_vect, gamma) <= 0) {
      0
    } else {
      uniroot(f = shifted_pearson_stat, interval = c(0, 0.5), n_vect, k_vect, pi_tilde_vect, gamma)$root
    }
  })
}


flag_outlier_dispersions <- function(pilot_rho_tilde_per_amplicon, outlier_mad_thresh) {
  positive_pilot_rho_tilde <- pilot_rho_tilde_per_amplicon[pilot_rho_tilde_per_amplicon > 0]
  if (length(positive_pilot_rho_tilde) >= 1L) {
    my_median <- median(positive_pilot_rho_tilde)
    my_mad <- stats::mad(positive_pilot_rho_tilde)
    outlier_thresh <- my_median + outlier_mad_thresh * my_mad
    ok_disp <- pilot_rho_tilde_per_amplicon <= outlier_thresh
  } else {
    ok_disp <- rep(TRUE, length(pilot_rho_tilde_per_amplicon))
    my_median <- 0
    my_mad <- 0
    outlier_thresh <- NA
  }
  return(list(ok_disp = ok_disp, median_pilot_rho = my_median,
              mad_pilot_rho = my_mad, outlier_thresh = outlier_thresh))
}


remove_outlier_dispersions_and_shrink_to_median <- function(pilot_rho_tilde_per_amplicon, outlier_mad_thresh) {
  my_median <- median(pilot_rho_tilde_per_amplicon)
  my_mad <- stats::mad(pilot_rho_tilde_per_amplicon)
  outlier_thresh <- my_median + outlier_mad_thresh * my_mad
  rho_tilde_shrunk <- median(pilot_rho_tilde_per_amplicon[pilot_rho_tilde_per_amplicon < outlier_thresh])
  pilot_rho_tilde_per_amplicon[pilot_rho_tilde_per_amplicon < outlier_thresh] <- rho_tilde_shrunk
  return(pilot_rho_tilde_per_amplicon)
}


run_fisher_exact_test <- function(data_list, nominal_fdr = 0.1, alternative = "greater") {
  n_amplicons <- ncol(data_list$n_mat_trt)
  p_vals <- sapply(X = seq_len(n_amplicons), FUN = function(i) {
    n_trt <- sum(data_list$n_mat_trt[,i])
    k_trt <- sum(data_list$k_mat_trt[,i])
    n_cntrl <- sum(data_list$n_mat_cntrl[,i])
    k_cntrl <- sum(data_list$k_mat_cntrl[,i])
    mat <- matrix(data = c(k_trt, n_trt - k_trt, k_cntrl, n_cntrl - k_cntrl), nrow = 2L)
    fit <- fisher.test(mat, alternative = alternative)
    fit$p.value
  })
  significant <- p.adjust(p = p_vals, method = "BH") <= nominal_fdr
  to_return <- data.frame(amplicon_id = data_list$amplicon_ids,
                          p_value = p_vals, significant = significant)
}
