#' Generate synthetic amplicon-seq data
#'
#' @param p number of amplicons
#' @param r number of treatment/control samples
#' @param pi_cntrl ([0,1] scalar) background editing rate to use across all control samples
#' @param editing_rate ([0,1] scalar) editing rate in samples where there is nonzero editing
#' @param n_amplicons_nonzero_editing (integer) number of amplicons
#' @param sample_size_mu (positive scalar) mean sample size
#' @param sample_size_theta (positive scalar) size (i.e., overdispersion) parameter for sample size
#' @param beta_binom_rho (positive scalar) dispersion parameter for the beta-binomial distribution
#'
#' @returns
#' @export
#'
#' @examples
generate_synthetic_amplicon_seq_data <- function(p, r, pi_cntrl, editing_rate, n_amplicons_nonzero_editing,
                                                 beta_binom_rho, sample_size_mu = 100000L, sample_size_theta = 15L) {
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
  generate_k_mat <- function(n_mat, pi_vect, beta_binom_rho, p, r) {
    k_mat <- sapply(X = seq_len(p), FUN = function(j) {
      VGAM::rbetabinom(n = r, size = n_mat[,j], prob = pi_vect[j], rho = beta_binom_rho)
    })
    return(k_mat)
  }
  k_mat_trt <- generate_k_mat(n_mat_trt, pi_trt_vect, beta_binom_rho, p, r)
  k_mat_cntrl <- generate_k_mat(n_mat_cntrl, pi_cntrl_vect, beta_binom_rho, p, r)

  # set column names
  colnames(n_mat_trt) <- colnames(n_mat_cntrl) <- colnames(k_mat_trt) <- colnames(k_mat_cntrl) <- seq(1L, p)
  out_list <- list(n_mat_trt = n_mat_trt,
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
#' @param data_list a list containing n_mat_trt, n_mat_cntrl, k_mat_trt, k_mat_cntrl
#' @param editing_threshold (scalar [0,1]) test whether editing is greater than this threshold
#'
#' @returns a data frame containing
#' @export
#'
#' @examples
#' data_list <- generate_synthetic_amplicon_seq_data(p = 20L, r = 3L, pi_cntrl = 0.05, editing_rate = 0.1,
#'                                                   n_amplicons_nonzero_editing = 2L, beta_binom_rho = 0.01)
#' result_df <- run_amplicon_seq_analysis(data_list) |>
#'   dplyr::mutate(amplicon_id = factor(amplicon_id, labels = seq_len(p), levels = seq_len(p)))
#' create_amplicon_seq_ci_plot(result_df)
run_amplicon_seq_analysis <- function(data_list, editing_threshold = 0.01, nominal_ci_coverage = 0.95, nominal_fdr = 0.1, rho = NULL, tail = "right") {
  r <- nrow(data_list$n_mat_trt)
  p <- ncol(data_list$n_mat_trt)

  # first, estimate of pi for each (amplicon, condition) pair
  estimate_pi_per_amplicon <- function(n_mat, k_mat) colSums(k_mat)/colSums(n_mat)
  pi_hat_trt_per_amplicon <- estimate_pi_per_amplicon(n_mat = data_list$n_mat_trt, k_mat = data_list$k_mat_trt)
  pi_hat_cntrl_per_amplicon <- estimate_pi_per_amplicon(n_mat = data_list$n_mat_cntrl, k_mat = data_list$k_mat_cntrl)

  # second, estimate rho via method of moments estimator
  if (is.null(rho)) {
    n_vect <- c(as.integer(data_list$n_mat_trt), as.integer(data_list$n_mat_cntrl))
    k_vect <- c(as.integer(data_list$k_mat_trt), as.integer(data_list$k_mat_cntrl))
    pi_hat_vect <- c(rep(pi_hat_trt_per_amplicon, each = r), rep(pi_hat_cntrl_per_amplicon, each = r))
    gamma <- 2 * p * (r - 1)
    shifted_pearson_stat <- function(cand_rho, n_vect, k_vect, pi_hat_vect, gamma) {
      sum((k_vect - n_vect * pi_hat_vect)^2/(n_vect * pi_hat_vect * (1 - pi_hat_vect) * (1 + (n_vect - 1) *  cand_rho))) - gamma
    }
    rho_hat <- uniroot(f = shifted_pearson_stat, interval = c(0, 0.5), n_vect, k_vect, pi_hat_vect, gamma)$root
  } else{
    rho_hat <- rho
  }

  # third, compute the standard error of the pi_hats
  compute_pi_hat_ses <- function(n_mat, pi_hat_per_amplicon, rho_hat) {
    sapply(X = seq_len(p), FUN = function(i) {
      n <- n_mat[,i]
      n_tot <- sum(n)
      pi_hat <- pi_hat_per_amplicon[[i]]
      sqrt(sum(n * pi_hat * (1 - pi_hat) * (1 + (n - 1) * rho_hat)))/n_tot
    })
  }
  pi_hat_se_trt_per_amplicon <- compute_pi_hat_ses(n_mat = data_list$n_mat_trt,
                                                   pi_hat_per_amplicon = pi_hat_trt_per_amplicon,
                                                   rho_hat = rho_hat)
  pi_hat_se_cntrl_per_amplicon <- compute_pi_hat_ses(n_mat = data_list$n_mat_cntrl,
                                                     pi_hat_per_amplicon = pi_hat_cntrl_per_amplicon,
                                                     rho_hat = rho_hat)

  # fourth, compute the theta_hats (editing rate estimates)
  theta_hat_per_amplicon <- (pi_hat_trt_per_amplicon - pi_hat_cntrl_per_amplicon)/(1 - pi_hat_cntrl_per_amplicon)

  # fifth, compute the theta_hat standard errors
  theta_hat_se_per_amplicon <- sapply(X = seq_len(p), FUN = function(i) {
    pi_hat_trt <- pi_hat_trt_per_amplicon[[i]]
    pi_hat_cntrl <- pi_hat_cntrl_per_amplicon[[i]]
    pi_hat_se_trt <- pi_hat_se_trt_per_amplicon[[i]]
    pi_hat_se_cntrl <- pi_hat_se_cntrl_per_amplicon[[i]]
    1/(1 - pi_hat_cntrl)^2 * sqrt((1 - pi_hat_cntrl)^2 * pi_hat_se_trt^2 + (1 - pi_hat_trt)^2 * pi_hat_se_cntrl^2)
  })

  # sixth, compute confidence intervals and standard errors
  mult_factor <- qnorm(p = (1 - nominal_ci_coverage)/2, lower.tail = FALSE)
  lower_theta_ci_per_amplicon <- pmax(pmin(theta_hat_per_amplicon - mult_factor * theta_hat_se_per_amplicon, 1), 0)
  upper_theta_ci_per_amplicon <- pmax(pmin(theta_hat_per_amplicon + mult_factor * theta_hat_se_per_amplicon, 1), 0)
  z_score_per_amplicon <- (theta_hat_per_amplicon - editing_threshold)/theta_hat_se_per_amplicon
  p_val_per_amplicon <- if (tail == "right") {
    pnorm(q = z_score_per_amplicon, lower.tail = FALSE)
  } else if (tail == "both") {
    2 * pmin(pnorm(q = z_score_per_amplicon, lower.tail = FALSE), pnorm(q = z_score_per_amplicon, lower.tail = TRUE))
  } else {
    stop("Tail not recognized.")
  }

  significant <- p.adjust(p = p_val_per_amplicon, method = "BH") <= nominal_fdr
  pi_hat_trt <- pi_hat_trt_per_amplicon
  pi_hat_cntrl <- pi_hat_cntrl_per_amplicon

  # return output
  to_return <- data.frame(amplicon_id = colnames(data_list$n_mat_trt),
                          theta_hat = theta_hat_per_amplicon,
                          theta_hat_clipped = pmax(pmin(theta_hat_per_amplicon, 1), 0),
                          theta_hat_se = theta_hat_se_per_amplicon,
                          theta_hat_lower_ci = lower_theta_ci_per_amplicon,
                          theta_hat_upper_ci = upper_theta_ci_per_amplicon,
                          pi_hat_trt = pi_hat_trt,
                          pi_hat_cntrl = pi_hat_cntrl,
                          pi_hat_se_trt = pi_hat_se_trt_per_amplicon,
                          pi_hat_se_cntrl = pi_hat_se_cntrl_per_amplicon,
                          p_value = pmax(p_val_per_amplicon, 1e-250),
                          significant = significant)
  rownames(to_return) <- NULL
  return(list(result_df = to_return, rho_hat = rho_hat))
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
  to_return <- data.frame(amplicon_id = colnames(data_list$n_mat_trt),
                          p_value = p_vals, significant = significant)
}

create_amplicon_seq_ci_plot <- function(result_df) {
  my_plot <- ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = amplicon_id, y = theta_hat_clipped)) +
    ggplot2::geom_point() + ggplot2::theme_bw() +
    ggplot2::geom_errorbar(mapping = ggplot2::aes(ymin = theta_hat_lower_ci,
                                                  ymax = theta_hat_upper_ci, width = 0)) +
    ggplot2::xlab("Amplicon") + ggplot2::ylab("Estimated editing rate")
}

create_amplicon_seq_p_value_plot <- function(result_df, min_p_value = 1e-250) {
  ggplot2::ggplot(data = result_df |>
                    dplyr::mutate(Significant = significant,
                                  p_value = pmax(p_value, min_p_value)),
         mapping = ggplot2::aes(x = amplicon_id, y = p_value, col = Significant)) +
    ggplot2::geom_point() + ggplot2::theme_bw(base_size = 9) + ggplot2::scale_color_manual(values = c("black", "dodgerblue2")) +
    ggplot2::xlab("Amplicon") + ggplot2::scale_y_continuous(trans = sceptre:::revlog_trans()) +
    ggplot2::ylab("p-value")
}

convert_data_list_into_mutation_frac_df <- function(data_list) {
  compute_mutation_frac <- function(n_mat, k_mat) {
    as.data.frame(k_mat/n_mat) |>
      tidyr::pivot_longer(cols = tidyr::everything(),
                          names_to = "amplicon", values_to = "mutation_frac") |>
      dplyr::arrange(amplicon)
  }
  mutation_frac_df <- rbind(compute_mutation_frac(data_list$n_mat_trt, data_list$k_mat_trt) |>
                              dplyr::mutate(Condition = "Treated"),
                            compute_mutation_frac(data_list$n_mat_cntrl, data_list$k_mat_cntrl) |>
                              dplyr::mutate(Condition = "Control"))

  n_vect <- colSums(data_list$n_mat_trt) + colSums(data_list$n_mat_cntrl)
  n_reads_df <- data.frame(amplicon = names(n_vect), n_reads = setNames(n_vect, NULL))

  return(list(mutation_frac_df = mutation_frac_df, n_reads_df = n_reads_df))
}

make_mutation_frac_plot <- function(to_plot) {
  ggplot2::ggplot(data = to_plot,
                  mapping = ggplot2::aes(x = Condition, y = mutation_frac, col = Condition)) +
    ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::xlab("Sample") + ggplot2::ylab("Mutation fraction") +
    ggplot2::facet_wrap(. ~ amplicon) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   strip.text = element_text(size = 8), legend.position = "bottom") +
    ggplot2::scale_color_manual(values = c("firebrick1", "dodgerblue1"))
}

make_n_reads_plot <- function(to_plot) {
  ggplot2::ggplot(data = to_plot,
                  mapping = ggplot2::aes(x = amplicon, y = n_reads)) +
    ggplot2::geom_point() + ggplot2::theme_bw()
}
