#' Generate synthetic amplicon-seq data
#'
#' @param p number of amplicons
#' @param r number of treatment/control samples
#' @param pi_cntrl ([0,1] scalar) background editing rate to use across all control samples
#' @param editing_rate ([0,1] scalar) editing rate in samples where there is nonzero editing
#' @param n_amplicons_nonzero_editing (integer) number of amplicons
#' @param amplicon_ids vector of amplicon identifiers
#' @param sample_size_mu (positive scalar) mean sample size
#' @param sample_size_theta (positive scalar) size (i.e., overdispersion) parameter for sample size
#' @param beta_binom_rho (positive scalar or vector) dispersion parameter for the beta-binomial distribution
#'
#' @returns A list containing the simulated treated/control count matrices, amplicon identifiers,
#' the generating mutation rates, and the beta-binomial dispersion values used in the simulation.
#' @export
#'
#' @examples
#' p <- 50L
#' sim_data <- generate_synthetic_amplicon_seq_data(
#'   p = 50L,
#'   r = 3L,
#'   pi_cntrl = 0.02,
#'   editing_rate = 0.15,
#'   n_amplicons_nonzero_editing = 2L,
#'   beta_binom_rho = 0.005,
#'   amplicon_ids = paste0("amplicon_", seq_len(p)),
#'   sample_size_mu = 50L
#' )
generate_synthetic_amplicon_seq_data <- function(p, r, pi_cntrl, editing_rate, n_amplicons_nonzero_editing,
                                                 beta_binom_rho, amplicon_ids, sample_size_mu = 100000L,
                                                 sample_size_theta = 15L) {
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


resolve_amplicon_ids <- function(data_list, n_amplicons) {
  amplicon_ids <- data_list$amplicon_ids
  if (is.null(amplicon_ids) || length(amplicon_ids) == 0L) {
    amplicon_ids <- colnames(data_list$n_mat_trt)
  }
  if (is.null(amplicon_ids) || length(amplicon_ids) == 0L) {
    amplicon_ids <- as.character(seq_len(n_amplicons))
  }
  if (length(amplicon_ids) != n_amplicons) {
    stop("Length of `amplicon_ids` must match the number of matrix columns.")
  }
  amplicon_ids
}


initialize_dispersion_diagnostics <- function() {
  list(ok_disp = logical(0),
       median_pilot_rho = NA_real_,
       mad_pilot_rho = NA_real_,
       outlier_thresh = NA_real_,
       shared_rho_hat = NA_real_)
}


#' Run amplicon-seq analysis
#'
#' @param data_list a list containing n_mat_trt, n_mat_cntrl, k_mat_trt, k_mat_cntrl
#' @param editing_threshold (scalar [0,1]) test whether editing is greater than this threshold
#' @param nominal_ci_coverage (scalar [0,1]) nominal coverage for Wald confidence intervals
#' @param nominal_fdr (scalar [0,1]) false discovery rate target for BH adjustment
#' @param global_rho logical; if TRUE, estimate one shared dispersion across QC-passing amplicons
#' @param rho optional scalar dispersion to use instead of estimating it
#' @param tail one of "right" or "both"
#' @param outlier_mad_thresh positive scalar determining the pilot-dispersion outlier threshold
#' @param min_mutation_count minimum total mutated reads across treated and control samples required for Wald inference
#'
#' @returns A list with two components:
#' \describe{
#'   \item{result_df}{A data frame containing point estimates for all amplicons and Wald inferential summaries
#'   for the amplicons that pass the mutation-count QC screen.}
#'   \item{dispersion_diagnostics}{A list containing dataset-level dispersion diagnostics for plotting.}
#' }
#' @export
#'
#' @examples
#' set.seed(1)
#' p <- 50L
#' amplicon_ids <- factor(x = paste0("amplicon_", seq_len(p)), levels = paste0("amplicon_", seq_len(p)))
#' beta_binom_rho <- c(0.005, rep(5e-4, times = p - 1L))
#' data_list <- generate_synthetic_amplicon_seq_data(p = p, r = 3L, pi_cntrl = 0.05, editing_rate = 0.15,
#'                                                   n_amplicons_nonzero_editing = 2L, beta_binom_rho,
#'                                                   amplicon_ids = amplicon_ids, sample_size_mu = 50L)
#' res <- run_amplicon_seq_analysis(data_list)
#'
#' # make plots
#' make_amplicon_seq_ci_plot(res$result_df) |> plot()
#' make_amplicon_seq_p_value_plot(res$result_df) |> plot()
#' make_pilot_dispersion_plot(res) |> plot()
run_amplicon_seq_analysis <- function(data_list, editing_threshold = 0.001, nominal_ci_coverage = 0.95,
                                      nominal_fdr = 0.1, global_rho = FALSE, rho = NULL, tail = "right",
                                      outlier_mad_thresh = 4, min_mutation_count = 20L) {
  n_mat_trt <- data_list$n_mat_trt
  n_mat_cntrl <- data_list$n_mat_cntrl
  k_mat_trt <- data_list$k_mat_trt
  k_mat_cntrl <- data_list$k_mat_cntrl
  p <- ncol(n_mat_trt)
  amplicon_ids <- resolve_amplicon_ids(data_list = data_list, n_amplicons = p)

  # 0. perform qc on total mutated reads, but keep all amplicons in the output
  total_mutated_reads <- colSums(k_mat_trt) + colSums(k_mat_cntrl)
  passes_mutation_count_qc_v <- total_mutated_reads >= min_mutation_count

  # 1. estimate of pi for each (amplicon, condition) pair
  estimate_pi_per_amplicon <- function(n_mat, k_mat) colSums(k_mat)/colSums(n_mat)
  estimate_jeffreys_pi_per_amplicon <- function(n_mat, k_mat) (colSums(k_mat) + 0.5)/(colSums(n_mat) + 1)
  pi_hat_trt_per_amplicon <- estimate_pi_per_amplicon(n_mat = n_mat_trt, k_mat = k_mat_trt)
  pi_hat_cntrl_per_amplicon <- estimate_pi_per_amplicon(n_mat = n_mat_cntrl, k_mat = k_mat_cntrl)
  pi_tilde_trt_per_amplicon <- estimate_jeffreys_pi_per_amplicon(n_mat = n_mat_trt, k_mat = k_mat_trt)
  pi_tilde_cntrl_per_amplicon <- estimate_jeffreys_pi_per_amplicon(n_mat = n_mat_cntrl, k_mat = k_mat_cntrl)

  # 2. estimate rho via method-of-moments on the QC-passing amplicons
  rho_hat_per_amplicon <- rep(NA_real_, p)
  pilot_rho_hat_per_amplicon <- rep(NA_real_, p)
  dispersion_outlier <- rep(NA, p)
  dispersion_diagnostics <- initialize_dispersion_diagnostics()

  if (any(passes_mutation_count_qc_v)) {
    n_mat_trt_qc_subset <- n_mat_trt[, passes_mutation_count_qc_v, drop = FALSE]
    n_mat_cntrl_qc_subset <- n_mat_cntrl[, passes_mutation_count_qc_v, drop = FALSE]
    k_mat_trt_qc_subset <- k_mat_trt[, passes_mutation_count_qc_v, drop = FALSE]
    k_mat_cntrl_qc_subset <- k_mat_cntrl[, passes_mutation_count_qc_v, drop = FALSE]
    pi_tilde_trt_qc_subset <- pi_tilde_trt_per_amplicon[passes_mutation_count_qc_v]
    pi_tilde_cntrl_qc_subset <- pi_tilde_cntrl_per_amplicon[passes_mutation_count_qc_v]

    if (is.null(rho)) {
      if (global_rho) {
        rho_hat_qc_subset_final <- estimate_global_rho(
          n_mat_trt = n_mat_trt_qc_subset,
          n_mat_cntrl = n_mat_cntrl_qc_subset,
          k_mat_trt = k_mat_trt_qc_subset,
          k_mat_cntrl = k_mat_cntrl_qc_subset,
          pi_hat_trt_per_amplicon = pi_tilde_trt_qc_subset,
          pi_hat_cntrl_per_amplicon = pi_tilde_cntrl_qc_subset
        )
      } else {
        rho_hat_qc_subset_final <- pilot_rho_hat_qc_subset <- estimate_rho_per_amplicon(
          n_mat_trt = n_mat_trt_qc_subset,
          n_mat_cntrl = n_mat_cntrl_qc_subset,
          k_mat_trt = k_mat_trt_qc_subset,
          k_mat_cntrl = k_mat_cntrl_qc_subset,
          pi_hat_trt_per_amplicon = pi_tilde_trt_qc_subset,
          pi_hat_cntrl_per_amplicon = pi_tilde_cntrl_qc_subset
        )
        pilot_rho_hat_per_amplicon[passes_mutation_count_qc_v] <- pilot_rho_hat_qc_subset
        dispersion_diagnostics <- flag_outlier_dispersions(pilot_rho_hats = pilot_rho_hat_qc_subset,
                                                           outlier_mad_thresh = outlier_mad_thresh)
        qc_subset_non_outlier_v <- dispersion_diagnostics$ok_disp
        dispersion_outlier[passes_mutation_count_qc_v] <- !qc_subset_non_outlier_v

        if (any(qc_subset_non_outlier_v)) {
          updated_rho_hat_qc_non_outlier_subset <- estimate_global_rho(
            n_mat_trt = n_mat_trt_qc_subset[, qc_subset_non_outlier_v, drop = FALSE],
            n_mat_cntrl = n_mat_cntrl_qc_subset[, qc_subset_non_outlier_v, drop = FALSE],
            k_mat_trt = k_mat_trt_qc_subset[, qc_subset_non_outlier_v, drop = FALSE],
            k_mat_cntrl = k_mat_cntrl_qc_subset[, qc_subset_non_outlier_v, drop = FALSE],
            pi_hat_trt_per_amplicon = pi_tilde_trt_qc_subset[qc_subset_non_outlier_v],
            pi_hat_cntrl_per_amplicon = pi_tilde_cntrl_qc_subset[qc_subset_non_outlier_v]
          )
          rho_hat_qc_subset_final[qc_subset_non_outlier_v] <- updated_rho_hat_qc_non_outlier_subset
          dispersion_diagnostics$shared_rho_hat <- updated_rho_hat_qc_non_outlier_subset[1]
        }
      }
    } else {
      rho_hat_qc_subset_final <- rep(rho, sum(passes_mutation_count_qc_v))
      dispersion_diagnostics$shared_rho_hat <- rho
    }

    rho_hat_per_amplicon[passes_mutation_count_qc_v] <- rho_hat_qc_subset_final
  }

  # 3. compute standard errors using Jeffreys-regularized proportions
  compute_pi_hat_ses <- function(n_mat, pi_for_var_per_amplicon, rho_hat_per_amplicon) {
    sapply(X = seq_len(ncol(n_mat)), FUN = function(i) {
      n <- n_mat[,i]
      n_tot <- sum(n)
      pi_for_var <- pi_for_var_per_amplicon[[i]]
      rho_hat <- rho_hat_per_amplicon[i]
      sqrt(sum(n * pi_for_var * (1 - pi_for_var) * (1 + (n - 1) * rho_hat)))/n_tot
    })
  }
  pi_hat_se_trt_per_amplicon <- rep(NA_real_, p)
  pi_hat_se_cntrl_per_amplicon <- rep(NA_real_, p)
  if (any(passes_mutation_count_qc_v)) {
    pi_hat_se_trt_per_amplicon[passes_mutation_count_qc_v] <- compute_pi_hat_ses(
      n_mat = n_mat_trt[, passes_mutation_count_qc_v, drop = FALSE],
      pi_for_var_per_amplicon = pi_tilde_trt_per_amplicon[passes_mutation_count_qc_v],
      rho_hat_per_amplicon = rho_hat_per_amplicon[passes_mutation_count_qc_v]
    )
    pi_hat_se_cntrl_per_amplicon[passes_mutation_count_qc_v] <- compute_pi_hat_ses(
      n_mat = n_mat_cntrl[, passes_mutation_count_qc_v, drop = FALSE],
      pi_for_var_per_amplicon = pi_tilde_cntrl_per_amplicon[passes_mutation_count_qc_v],
      rho_hat_per_amplicon = rho_hat_per_amplicon[passes_mutation_count_qc_v]
    )
  }

  # 4. compute the theta_hats (editing rate estimates) for all amplicons
  theta_hat_per_amplicon <- (pi_hat_trt_per_amplicon - pi_hat_cntrl_per_amplicon)/(1 - pi_hat_cntrl_per_amplicon)

  # 5. compute the theta_hat standard errors and Wald inference for the QC-passing amplicons
  theta_hat_se_per_amplicon <- rep(NA_real_, p)
  lower_theta_ci_per_amplicon <- rep(NA_real_, p)
  upper_theta_ci_per_amplicon <- rep(NA_real_, p)
  p_val_per_amplicon <- rep(NA_real_, p)
  significant <- rep(FALSE, p)

  if (any(passes_mutation_count_qc_v)) {
    theta_hat_se_qc_subset <- sapply(X = which(passes_mutation_count_qc_v), FUN = function(i) {
      pi_tilde_trt <- pi_tilde_trt_per_amplicon[[i]]
      pi_tilde_cntrl <- pi_tilde_cntrl_per_amplicon[[i]]
      pi_hat_se_trt <- pi_hat_se_trt_per_amplicon[[i]]
      pi_hat_se_cntrl <- pi_hat_se_cntrl_per_amplicon[[i]]
      1/(1 - pi_tilde_cntrl)^2 *
        sqrt((1 - pi_tilde_cntrl)^2 * pi_hat_se_trt^2 + (1 - pi_tilde_trt)^2 * pi_hat_se_cntrl^2)
    })
    theta_hat_se_per_amplicon[passes_mutation_count_qc_v] <- theta_hat_se_qc_subset

    mult_factor <- qnorm(p = (1 - nominal_ci_coverage)/2, lower.tail = FALSE)
    lower_theta_ci_per_amplicon[passes_mutation_count_qc_v] <-
      pmax(pmin(theta_hat_per_amplicon[passes_mutation_count_qc_v] - mult_factor * theta_hat_se_qc_subset, 1), 0)
    upper_theta_ci_per_amplicon[passes_mutation_count_qc_v] <-
      pmax(pmin(theta_hat_per_amplicon[passes_mutation_count_qc_v] + mult_factor * theta_hat_se_qc_subset, 1), 0)
    z_score_qc_subset <- (theta_hat_per_amplicon[passes_mutation_count_qc_v] - editing_threshold)/theta_hat_se_qc_subset
    p_val_qc_subset <- if (tail == "right") {
      pnorm(q = z_score_qc_subset, lower.tail = FALSE)
    } else if (tail == "both") {
      2 * pmin(pnorm(q = z_score_qc_subset, lower.tail = FALSE),
               pnorm(q = z_score_qc_subset, lower.tail = TRUE))
    } else {
      stop("Tail not recognized.")
    }
    p_val_per_amplicon[passes_mutation_count_qc_v] <- pmax(p_val_qc_subset, 1e-250)
    significant[passes_mutation_count_qc_v] <-
      p.adjust(p = p_val_per_amplicon[passes_mutation_count_qc_v], method = "BH") <= nominal_fdr
  }

  # return output
  to_return <- data.frame(amplicon_id = amplicon_ids,
                          total_mutated_reads = total_mutated_reads,
                          passes_mutation_count_qc = passes_mutation_count_qc_v,
                          theta_hat = theta_hat_per_amplicon,
                          rho_hat = rho_hat_per_amplicon,
                          pilot_rho_hat = pilot_rho_hat_per_amplicon,
                          dispersion_outlier = dispersion_outlier,
                          theta_hat_clipped = pmax(pmin(theta_hat_per_amplicon, 1), 0),
                          theta_hat_se = theta_hat_se_per_amplicon,
                          theta_hat_lower_ci = lower_theta_ci_per_amplicon,
                          theta_hat_upper_ci = upper_theta_ci_per_amplicon,
                          pi_hat_trt = pi_hat_trt_per_amplicon,
                          pi_hat_cntrl = pi_hat_cntrl_per_amplicon,
                          pi_hat_se_trt = pi_hat_se_trt_per_amplicon,
                          pi_hat_se_cntrl = pi_hat_se_cntrl_per_amplicon,
                          p_value = p_val_per_amplicon,
                          significant = significant)
  rownames(to_return) <- NULL
  out <- list(result_df = to_return,
              dispersion_diagnostics = dispersion_diagnostics)
  out
}


estimate_global_rho <- function(n_mat_trt, n_mat_cntrl, k_mat_trt, k_mat_cntrl,
                                pi_hat_trt_per_amplicon, pi_hat_cntrl_per_amplicon) {
  n_vect <- c(as.integer(n_mat_trt), as.integer(n_mat_cntrl))
  k_vect <- c(as.integer(k_mat_trt), as.integer(k_mat_cntrl))
  p <- ncol(n_mat_trt)
  r <- nrow(n_mat_trt)
  pi_hat_vect <- c(rep(pi_hat_trt_per_amplicon, each = r), rep(pi_hat_cntrl_per_amplicon, each = r))
  gamma <- 2 * p * (r - 1)
  shifted_pearson_stat <- function(cand_rho, n_vect, k_vect, pi_hat_vect, gamma) {
    sum((k_vect - n_vect * pi_hat_vect)^2/(n_vect * pi_hat_vect * (1 - pi_hat_vect) * (1 + (n_vect - 1) *  cand_rho))) - gamma
  }
  rho_hat <- if (shifted_pearson_stat(0, n_vect, k_vect, pi_hat_vect, gamma) <= 0) {
    0
  } else {
    uniroot(f = shifted_pearson_stat, interval = c(0, 0.5), n_vect, k_vect, pi_hat_vect, gamma)$root
  }
  rho_hat_per_amplicon <- rep(rho_hat, p)
  return(rho_hat_per_amplicon)
}


estimate_rho_per_amplicon <- function(n_mat_trt, n_mat_cntrl, k_mat_trt, k_mat_cntrl,
                                      pi_hat_trt_per_amplicon, pi_hat_cntrl_per_amplicon) {
  p <- ncol(n_mat_trt)
  r <- nrow(n_mat_trt)
  sapply(X = seq_len(p), FUN = function(amplicon_idx) {
    n_vect <- c(n_mat_trt[,amplicon_idx], n_mat_cntrl[,amplicon_idx])
    k_vect <- c(k_mat_trt[,amplicon_idx], k_mat_cntrl[,amplicon_idx])
    gamma <- 2 * r - 2
    pi_hat_vect <- c(rep(pi_hat_trt_per_amplicon[amplicon_idx], each = r),
                     rep(pi_hat_cntrl_per_amplicon[amplicon_idx], each = r))
    shifted_pearson_stat <- function(cand_rho, n_vect, k_vect, pi_hat_vect, gamma) {
      sum( (k_vect - n_vect * pi_hat_vect)^2/(n_vect * pi_hat_vect * (1 - pi_hat_vect) * (1 + (n_vect - 1) *  cand_rho)) ) - gamma
    }
    if (shifted_pearson_stat(0, n_vect, k_vect, pi_hat_vect, gamma) <= 0) {
      0
    } else {
      uniroot(f = shifted_pearson_stat, interval = c(0, 0.5), n_vect, k_vect, pi_hat_vect, gamma)$root
    }
  })
}


flag_outlier_dispersions <- function(pilot_rho_hats, outlier_mad_thresh) {
  my_median <- median(pilot_rho_hats)
  my_mad <- stats::mad(pilot_rho_hats)
  outlier_thresh <- my_median + outlier_mad_thresh * my_mad
  ok_disp <- pilot_rho_hats <= outlier_thresh
  return(list(ok_disp = ok_disp, median_pilot_rho = my_median,
              mad_pilot_rho = my_mad, outlier_thresh = outlier_thresh))
}


remove_outlier_dispersions_and_shrink_to_median <- function(pilot_rho_hat_per_amplicon, outlier_mad_thresh) {
  my_median <- median(pilot_rho_hat_per_amplicon)
  my_mad <- stats::mad(pilot_rho_hat_per_amplicon)
  outlier_thresh <- my_median + outlier_mad_thresh * my_mad
  rho_hat_shrunk <- median(pilot_rho_hat_per_amplicon[pilot_rho_hat_per_amplicon < outlier_thresh])
  pilot_rho_hat_per_amplicon[pilot_rho_hat_per_amplicon < outlier_thresh] <- rho_hat_shrunk
  return(pilot_rho_hat_per_amplicon)
}


run_fisher_exact_test <- function(data_list, nominal_fdr = 0.1, alternative = "greater") {
  n_amplicons <- ncol(data_list$n_mat_trt)
  amplicon_ids <- resolve_amplicon_ids(data_list = data_list, n_amplicons = n_amplicons)
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
  to_return <- data.frame(amplicon_id = amplicon_ids,
                          p_value = p_vals, significant = significant)
}
