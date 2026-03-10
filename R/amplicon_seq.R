#' Generate synthetic amplicon-seq data
#'
#' @param p number of amplicons
#' @param r number of treatment/control samples
#' @param pi_cntrl ([0,1] scalar) background editing rate to use across all control samples
#' @param editing_rate ([0,1] scalar) editing rate in samples where there is nonzero editing
#' @param n_amplicons_nonzero_editing (integer) number of amplicons
#' @param sample_size_mu (positive scalar) mean sample size
#' @param sample_size_theta (positive scalar) size (i.e., overdispersion) parameter for sample size
#' @param beta_binom_rho (positive scalar or vector) dispersion parameter for the beta-binomial distribution
#'
#' @returns
#' @export
#'
#' @examples
generate_synthetic_amplicon_seq_data <- function(p, r, pi_cntrl, editing_rate, n_amplicons_nonzero_editing,
                                                 beta_binom_rho, amplicon_ids, sample_size_mu = 100000L, sample_size_theta = 15L) {
  pi_trt_under_editing <- (pi_cntrl + editing_rate - pi_cntrl * editing_rate)
  pi_cntrl_vect <- rep(pi_cntrl, p)
  pi_trt_vect <- c(rep(pi_trt_under_editing, n_amplicons_nonzero_editing), rep(pi_cntrl, p - n_amplicons_nonzero_editing))
  df_list <- list()
  counter <- 1L
  n_trt_samples <- p * r

  # generate n_mats
  n_mat_list <- lapply(X = seq_len(2L), FUN = function(i) {
    matrix(data = MASS::rnegbin(n = n_trt_samples,
                                mu = sample_size_mu,
                                theta = sample_size_theta), nrow = r, ncol = p)
  })
  n_mat_trt <- n_mat_list[[1]]; n_mat_cntrl <- n_mat_list[[2]]

  # generate k_mats
  if (length(beta_binom_rho) == 1L) beta_binom_rho <- rep(beta_binom_rho, p)
  generate_k_mat <- function(n_mat, pi_vect, beta_binom_rho, p, r) {
    k_mat <- sapply(X = seq_len(p), FUN = function(j) {
      VGAM::rbetabinom(n = r, size = n_mat[,j], prob = pi_vect[j], rho = beta_binom_rho[j])
    })
    return(k_mat)
  }
  k_mat_trt <- generate_k_mat(n_mat_trt, pi_trt_vect, beta_binom_rho, p, r)
  k_mat_cntrl <- generate_k_mat(n_mat_cntrl, pi_cntrl_vect, beta_binom_rho, p, r)

  # set column names
  out_list <- list(amplicon_ids = amplicon_ids,
                   n_mat_trt = n_mat_trt,
                   n_mat_cntrl = n_mat_cntrl,
                   k_mat_trt = k_mat_trt,
                   k_mat_cntrl = k_mat_cntrl,
                   pi_trt_vect = pi_trt_vect,
                   pi_cntrl_vect = pi_cntrl_vect,
                   beta_binom_rho = beta_binom_rho)
  return(out_list)
}


#' Run amplicon-seq analysis
#'
#' @param data_list a list containing n_mat_trt, n_mat_cntrl, k_mat_trt, k_mat_cntrl, and amplicon_ids
#' @param editing_threshold (scalar [0,1]) test whether editing is greater than this threshold
#'
#' @returns a data frame containing
#' @export
#'
#' @examples
#' set.seed(1)
#' p <- 25L
#' amplicon_ids <- factor(x = paste0("amplicon_", seq_len(p)), levels = paste0("amplicon_", seq_len(p)))
#' beta_binom_rho <- c(0.005, rep(5e-4, times = p - 1L))
#' data_list <- generate_synthetic_amplicon_seq_data(p = p, r = 3L, pi_cntrl = 0.05, editing_rate = 0.15,
#'                                                   n_amplicons_nonzero_editing = 2L, beta_binom_rho,
#'                                                   amplicon_ids = amplicon_ids)
#' res <- run_amplicon_seq_analysis(data_list)
#'
#' # make plots
#' make_amplicon_seq_ci_plot(res$result_df) |> plot()
#' make_amplicon_seq_p_value_plot(res$result_df) |> plot()
#' make_pilot_dispersion_plot(res) |> plot()
run_amplicon_seq_analysis <- function(data_list, editing_threshold = 0.001, nominal_ci_coverage = 0.95,
                                      nominal_fdr = 0.1, global_rho = FALSE, rho = NULL, tail = "right",
                                      outlier_mad_thresh = 4, min_mutation_count = 10L) {
  # 0. extract the data
  n_mat_trt <- data_list$n_mat_trt
  n_mat_cntrl <- data_list$n_mat_cntrl
  k_mat_trt <- data_list$k_mat_trt
  k_mat_cntrl <- data_list$k_mat_cntrl

  # 1. estimate of pi and theta
  estimate_pi_per_amplicon <- function(n_mat, k_mat) colSums(k_mat)/colSums(n_mat)
  pi_hat_trt_per_amplicon <- estimate_pi_per_amplicon(n_mat = n_mat_trt, k_mat = k_mat_trt)
  pi_hat_cntrl_per_amplicon <- estimate_pi_per_amplicon(n_mat = n_mat_cntrl, k_mat = k_mat_cntrl)
  theta_hat_per_amplicon <- (pi_hat_trt_per_amplicon - pi_hat_cntrl_per_amplicon)/(1 - pi_hat_cntrl_per_amplicon)

  # 2. perform qc, retaining only the amplicons with a mutation count of at least `min_mutation_count` across samples
  mut_count_ok_v <- colSums(k_mat_trt) + colSums(k_mat_cntrl) >= min_mutation_count
  if (any(!mut_count_ok_v)) {
    # construct fail qc result data frame
    to_return_fail_qc <- data.frame(amplicon_id = data_list$amplicon_ids[!mut_count_ok_v],
                                    theta_hat = theta_hat_per_amplicon[!mut_count_ok_v],
                                    rho_hat = NA, pilot_rho_hat = NA,
                                    dispersion_outlier = NA,
                                    theta_hat_clipped = pmax(pmin(theta_hat_per_amplicon[!mut_count_ok_v], 1), 0),
                                    theta_hat_se = NA, theta_hat_lower_ci = NA,
                                    theta_hat_upper_ci = NA,
                                    pi_hat_trt = pi_hat_trt_per_amplicon[!mut_count_ok_v],
                                    pi_hat_cntrl = pi_hat_cntrl_per_amplicon[!mut_count_ok_v],
                                    pi_hat_se_trt = NA,
                                    pi_hat_se_cntrl = NA, p_value = NA, significant = NA, pass_qc = FALSE)
    # if all fail qc, return this data frame
    if (all(!mut_count_ok_v)) {
      ret <- list(result_df = to_return_fail_qc)
      return(ret)
    }
  }

  # 3. compute pi tilde, or the Jeffrey-regularized pi estimator
  n_mat_trt <- n_mat_trt[, mut_count_ok_v, drop = FALSE]
  n_mat_cntrl <- n_mat_cntrl[, mut_count_ok_v, drop = FALSE]
  k_mat_trt <- k_mat_trt[, mut_count_ok_v, drop = FALSE]
  k_mat_cntrl <- k_mat_cntrl[, mut_count_ok_v, drop = FALSE]
  estimate_reg_pi_per_amplicon <- function(n_mat, k_mat) (colSums(k_mat) + 0.5)/(colSums(n_mat) + 1)
  pi_tilde_trt_per_amplicon <- estimate_reg_pi_per_amplicon(n_mat = n_mat_trt, k_mat = k_mat_trt)
  pi_tilde_cntrl_per_amplicon <- estimate_reg_pi_per_amplicon(n_mat = n_mat_cntrl, k_mat = k_mat_cntrl)

  # 4. estimate rho via method of moments estimator
  pilot_rho_tilde_per_amplicon <- NA
  dispersion_outlier <- NA
  dispersion_diagnostics <- NULL
  if (is.null(rho)) {
    if (global_rho) {
      rho_tilde_per_amplicon <- estimate_global_rho(n_mat_trt = n_mat_trt, n_mat_cntrl = n_mat_cntrl,
                                                    k_mat_trt = k_mat_trt, k_mat_cntrl = k_mat_cntrl,
                                                    pi_tilde_trt_per_amplicon = pi_tilde_trt_per_amplicon,
                                                    pi_tilde_cntrl_per_amplicon = pi_tilde_cntrl_per_amplicon)
    } else {
      # obtain pilot estimate of rho for each amplicon
      rho_tilde_per_amplicon <-
        pilot_rho_tilde_per_amplicon <- estimate_rho_per_amplicon(n_mat_trt = n_mat_trt, n_mat_cntrl = n_mat_cntrl,
                                                                  k_mat_trt = k_mat_trt, k_mat_cntrl = k_mat_cntrl,
                                                                  pi_tilde_trt_per_amplicon = pi_tilde_trt_per_amplicon,
                                                                  pi_tilde_cntrl_per_amplicon = pi_tilde_cntrl_per_amplicon)
      dispersion_diagnostics <- flag_outlier_dispersions(pilot_rho_tilde_per_amplicon = pilot_rho_tilde_per_amplicon,
                                                         outlier_mad_thresh = outlier_mad_thresh)
      disp_ok_v <- dispersion_diagnostics$ok_disp
      updated_rho_tilde_per_amplicon <- estimate_global_rho(n_mat_trt = n_mat_trt[,disp_ok_v, drop = FALSE],
                                                            n_mat_cntrl = n_mat_cntrl[,disp_ok_v, drop = FALSE],
                                                            k_mat_trt = k_mat_trt[,disp_ok_v, drop = FALSE],
                                                            k_mat_cntrl = k_mat_cntrl[,disp_ok_v, drop = FALSE],
                                                            pi_tilde_trt_per_amplicon = pi_tilde_trt_per_amplicon[disp_ok_v, drop = FALSE],
                                                            pi_tilde_cntrl_per_amplicon = pi_tilde_cntrl_per_amplicon[disp_ok_v, drop = FALSE])
      dispersion_diagnostics$shared_rho_hat <- updated_rho_tilde_per_amplicon[1]
      rho_tilde_per_amplicon[disp_ok_v] <- updated_rho_tilde_per_amplicon
      dispersion_outlier <- !disp_ok_v
    }
  } else {
    rho_tilde_per_amplicon <- rep(rho, ncol(n_mat_trt))
  }

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
  mult_factor <- qnorm(p = (1 - nominal_ci_coverage)/2, lower.tail = FALSE)
  lower_theta_ci_per_amplicon <- pmax(pmin(theta_hat_per_amplicon[mut_count_ok_v] - mult_factor * theta_tilde_se_per_amplicon, 1), 0)
  upper_theta_ci_per_amplicon <- pmax(pmin(theta_hat_per_amplicon[mut_count_ok_v] + mult_factor * theta_tilde_se_per_amplicon, 1), 0)
  z_score_per_amplicon <- (theta_hat_per_amplicon[mut_count_ok_v] - editing_threshold)/theta_tilde_se_per_amplicon
  p_val_per_amplicon <- if (tail == "right") {
    pnorm(q = z_score_per_amplicon, lower.tail = FALSE)
  } else if (tail == "both") {
    2 * pmin(pnorm(q = z_score_per_amplicon, lower.tail = FALSE), pnorm(q = z_score_per_amplicon, lower.tail = TRUE))
  } else {
    stop("Tail not recognized.")
  }
  significant <- p.adjust(p = p_val_per_amplicon, method = "BH") <= nominal_fdr

  # 8. return output
  to_return <- data.frame(amplicon_id = data_list$amplicon_ids[mut_count_ok_v],
                          theta_hat = theta_hat_per_amplicon[mut_count_ok_v],
                          rho_hat = rho_tilde_per_amplicon,
                          pilot_rho_hat = pilot_rho_tilde_per_amplicon,
                          dispersion_outlier = dispersion_outlier,
                          theta_hat_clipped = pmax(pmin(theta_hat_per_amplicon[mut_count_ok_v], 1), 0),
                          theta_hat_se = theta_tilde_se_per_amplicon,
                          theta_hat_lower_ci = lower_theta_ci_per_amplicon,
                          theta_hat_upper_ci = upper_theta_ci_per_amplicon,
                          pi_hat_trt = pi_hat_trt_per_amplicon[mut_count_ok_v],
                          pi_hat_cntrl = pi_hat_cntrl_per_amplicon[mut_count_ok_v],
                          pi_hat_se_trt = pi_tilde_se_trt_per_amplicon,
                          pi_hat_se_cntrl = pi_tilde_se_cntrl_per_amplicon,
                          p_value = pmax(p_val_per_amplicon, 1e-250),
                          significant = significant,
                          pass_qc = TRUE)
    if (any(!mut_count_ok_v)) to_return <- rbind(to_return, to_return_fail_qc)
    out <- list(result_df = to_return,
                dispersion_diagnostics = dispersion_diagnostics[-1])
    return(out)
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
  rho_tilde_per_amplicon <- rep(rho_tilde, p)
  return(rho_tilde_per_amplicon)
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
