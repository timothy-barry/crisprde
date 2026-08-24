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


fit_multirep_guideseq_occupancy <- function(Y_mat, incorporate_occupancy_info = TRUE) {
  MIN_NONZERO_COUNT <- 25L
  if (is.null(colnames(Y_mat))) warning("Y_mat must have column names (to identify the windows).")
  X <- Y_mat > 0
  storage.mode(X) <- "integer"

  nonzero_replicate_count <- rowSums(X)
  if (any(nonzero_replicate_count <= MIN_NONZERO_COUNT)) {
    offending_rows <- paste0(which(nonzero_replicate_count <= MIN_NONZERO_COUNT), collapse = ", ")
    msg <- paste0("Row ",  offending_rows, " has fewer than ", MIN_NONZERO_COUNT, " windows with a nonzero count. Consider dropping this sample or combining this sample with another (e.g., by pooling together primer chanels within a replicate).")
    warning(msg)
  }

  Omega <- as.matrix(generate_omega(nrow(Y_mat)))
  col_keys <- apply(X, 2, paste0, collapse = "")
  pi_hat <- NULL
  tbp_pattern_df <- NULL
  occupancy_pattern_map <- NULL

  if (incorporate_occupancy_info) {
    pi_hat <- fit_tbp_model(X)
    if (!is.null(pi_hat)) {
      tbp_pattern_df <- pmf_tbp(pi_hat)
    } else {
      warning("Cannot fit occupancy model; defaulting to count-only model.")
      incorporate_occupancy_info <- FALSE
    }
  }
  omega_keys <- apply(Omega, 1, paste0, collapse = "")
  occupancy_pattern_map <- match(col_keys, omega_keys)

  ret <- list(X = X,
              Omega = Omega,
              col_keys = col_keys,
              incorporate_occupancy_info = incorporate_occupancy_info,
              pi_hat = pi_hat,
              tbp_pattern_df = tbp_pattern_df,
              occupancy_pattern_map = occupancy_pattern_map)
  return(ret)
}

fit_multirep_guideseq_count_null <- function(Y_mat, c_tukey_beta = 5, c_tukey_sigma = 5, robust_fit = TRUE) {
  mu_theta_hat_mat <- apply(X = Y_mat, MARGIN = 1, FUN = function(curr_row) {
    y_plus <- curr_row[curr_row > 0]
    if (robust_fit) {
      fit <- fit_rob_nb_univariate(y = y_plus - 1, c.tukey.beta = c_tukey_beta, c.tukey.sigma = c_tukey_sigma)
    } else {
      fit <- fit_nb_univariate(y = y_plus - 1)
    }
    fit[c("mu", "theta")]
  }) |> t()

  return(mu_theta_hat_mat)
}


get_p_values_given_test_stats_prob_vector <- function(test_stat_v_in, right_tail_prob_v_in) {
  p_vals <- numeric(length(test_stat_v_in))
  neg_idx <- which(test_stat_v_in < 0)
  ok_idx <- which(test_stat_v_in >= 0)
  p_vals[neg_idx] <- 1
  tail_idx <- as.integer(test_stat_v_in[ok_idx]) + 1L
  too_large <- tail_idx > length(right_tail_prob_v_in)
  p_vals[ok_idx[too_large]] <- right_tail_prob_v_in[length(right_tail_prob_v_in)]
  p_vals[ok_idx[!too_large]] <- right_tail_prob_v_in[tail_idx[!too_large]]
  return(p_vals)
}


score_multirep_guideseq_fit <- function(Y_mat, occupancy_fit, mu_theta_hat_mat, lambda = 20,
                                        multiplicity_alpha = 0.1, annotated_clustered_count_df = NULL) {
  # unpack items
  X <- occupancy_fit$X
  Omega <- occupancy_fit$Omega
  col_keys <- occupancy_fit$col_keys
  incorporate_occupancy_info <- occupancy_fit$incorporate_occupancy_info
  occupancy_pattern_map <- occupancy_fit$occupancy_pattern_map

  # compute total UMI count and initialize p-value vector
  total_umi_counts <- colSums(Y_mat)
  occupancy_counts <- colSums(X)

  # iterate over occupancy patterns
  if (!incorporate_occupancy_info) {
    # no occupancy info
    p_vals <- numeric(length = length(total_umi_counts))
    test_stats <- total_umi_counts - occupancy_counts
    right_tail_prob_list <- get_right_tail_prob_list(mu_theta_hat_mat = mu_theta_hat_mat,
                                                     max_needed = max(test_stats),
                                                     Omega = Omega)
    for (i in seq_along(right_tail_prob_list)) {
      idxs <- which(occupancy_pattern_map == i)
      p_vals[idxs] <- get_p_values_given_test_stats_prob_vector(
        test_stat_v_in = test_stats[idxs],
        right_tail_prob_v_in = right_tail_prob_list[[i]]
      )
    }
  } else {
    # with occupancy info -- mixture over the occupancy patterns
    log_pi_hat <- log(occupancy_fit$pi_hat)
    window_log_pi_sum <- as.numeric(crossprod(log_pi_hat, X))
    pattern_log_pi_sum <- as.numeric(Omega %*% log_pi_hat)
    test_stats <- (total_umi_counts - occupancy_counts) - lambda * window_log_pi_sum
    max_needed <- max(0L, ceiling(max(test_stats) + lambda * max(pattern_log_pi_sum)))
    right_tail_prob_list <- get_right_tail_prob_list(mu_theta_hat_mat = mu_theta_hat_mat,
                                                     max_needed = max_needed, Omega = Omega)
    l <- sapply(X = seq_len(nrow(Omega)), FUN = function(i) {
      sum_start <- ceiling(test_stats + lambda * pattern_log_pi_sum[i])
      nb_piece <- get_p_values_given_test_stats_prob_vector(
        test_stat_v_in = sum_start,
        right_tail_prob_v_in = right_tail_prob_list[[i]]
      )
      occupancy_fit$tbp_pattern_df$pmf[i] * nb_piece
    }, simplify = FALSE)
    p_vals <- Reduce(f = "+", x = l)
  }

  ##### SECOND PART OF FUNCTION
  p_vals <- pmin(1, p_vals)
  q_vals <- p.adjust(p = p_vals, method = "BH")
  nominated_window <- (q_vals < multiplicity_alpha)
  lambda_out <- if (incorporate_occupancy_info) lambda else NA

  res_df <- data.frame(window = colnames(Y_mat),
                       p_value = p_vals,
                       test_stat = test_stats,
                       nominated_window = nominated_window,
                       umi_count = total_umi_counts,
                       lambda = lambda_out,
                       occupancy_pattern = col_keys) |> dplyr::arrange(p_value)
  rownames(res_df) <- NULL
  if (!is.null(annotated_clustered_count_df)) {
    right_df <- annotated_clustered_count_df |>
      dplyr::select(window, dplyr::starts_with(c("homology", "window", "overlaps"))) |>
      dplyr::filter(window %in% res_df$window) |>
      dplyr::distinct()
    res_df <- dplyr::left_join(res_df, right_df, by = "window")
  }
  ests_list <- list(mu_theta_hat_mat = mu_theta_hat_mat)
  if (incorporate_occupancy_info) ests_list$pi_hat <- occupancy_fit$pi_hat
  ret <- list(res_df = res_df, ests_list = ests_list)
  return(ret)
}


get_right_tail_prob_list <- function(mu_theta_hat_mat, max_needed, Omega) {
  pmf_list <- apply(X = mu_theta_hat_mat, MARGIN = 1, FUN = function(curr_row) {
    max_count <- max(max_needed, qnbinom(p = 1e-16, mu = curr_row[["mu"]], size = curr_row[["theta"]], lower.tail = FALSE))
    dnbinom(x = seq(0L, max_count), mu = curr_row[["mu"]], size = curr_row[["theta"]])
  }, simplify = FALSE)
  conv_pmf_list <- apply(X = Omega, MARGIN = 1, FUN = function(curr_row) {
    convolve_pmf_list(pmf_list[as.logical(curr_row)])
  }, simplify = FALSE)
  right_tail_prob_list <- lapply(X = conv_pmf_list, FUN = function(curr_pmf) {
    rev(cumsum(rev(curr_pmf)))
  })
  return(right_tail_prob_list)
}



#' Run multivariate guide-seq method
#'
#' Runs the multivariate guide-seq method on a replicate-by-window matrix of counts
#'
#' @param Y_mat the r x m integer matrix of UMI counts, where r is the number of replicates and m is the number of windows
#' @param incorporate_occupancy_info a boolean (T/F) indicating whether to incorporate occupancy information into the p-value calculation
#' @param annotated_clustered_count_df optional output of `annotate_clustered_count_df()`; if supplied, homology annotations are joined to the result data frame by window
#'
#' @export
#' @returns a data frame containing a p-value for each window
#' @examples
#' set.seed(42)
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
#' m_alt <- 15
#' alt_dat <- simulate_multirep_guideseq_data(pi, mu_vect, theta_vect, m_alt)
#'
#' # COMBINED DATA
#' Y_mat <- cbind(alt_dat, null_dat)
#' colnames(Y_mat) <- paste0("window_", seq_len(ncol(Y_mat)))
#' incorporate_occupancy_info <- TRUE
#' multiplicity_alpha <- 0.2
#'
#' # RUN METHOD
#' res_df <- run_multireplicate_guideseq_method(Y_mat, incorporate_occupancy_info = TRUE)
#'
#' # EVALUATE RESULT
#' n_correct_nominations <- sum(res_df$nominated_window[seq(1, m_alt)])
#' n_total_nominations <- sum(res_df$nominated_window)
#' fdp <- (n_total_nominations - n_correct_nominations)/n_total_nominations
run_multireplicate_guideseq_method <- function(Y_mat, incorporate_occupancy_info = TRUE,
                                               multiplicity_alpha = 0.1, lambda = 20,
                                               c_tukey_beta = 5, c_tukey_sigma = 5,
                                               robust_fit = TRUE, annotated_clustered_count_df = NULL) {
  occupancy_fit <- fit_multirep_guideseq_occupancy(Y_mat = Y_mat,
                                                   incorporate_occupancy_info = incorporate_occupancy_info)
  mu_theta_hat_mat <- fit_multirep_guideseq_count_null(Y_mat = Y_mat,
                                                       c_tukey_beta = c_tukey_beta,
                                                       c_tukey_sigma = c_tukey_sigma,
                                                       robust_fit = robust_fit)
  ret <- score_multirep_guideseq_fit(Y_mat = Y_mat,
                                     occupancy_fit = occupancy_fit,
                                     mu_theta_hat_mat = mu_theta_hat_mat,
                                     lambda = lambda,
                                     multiplicity_alpha = multiplicity_alpha,
                                     annotated_clustered_count_df = annotated_clustered_count_df)
  return(ret)
}


#' Cluster loci
#'
#' Clusters loci via single-linkage clustering
#'
#' @param count_df a data frame containing columns chr, coord
#' @param thresh clustering threshold; occupied bases within this distance are clustered together
#' @param padding amount of padding to add to either side of each cluster
#'
#' @returns `count_df` with the additional columns appended:
#' - `window` (a string indicating the cluster to which a given base belongs)
#' - `cluster_chr` (chromosome of the group)
#' - `min_cluster_coord` (minimum coordinate of the group)
#' - `max_cluster_coord` (maximum coordinate of the group)
#' @export
#'
#' @examples
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' count_df <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#'  dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & treated & replicate_id %in% 1:2) |>
#'  dplyr::filter(chr != "chrM") |>
#'  dplyr::select(chr, coord, strand, umi_count, primer_type, replicate_id)
#' clustered_count_df <- cluster_loci(count_df)
cluster_loci <- function(count_df, thresh = 100L, padding = 5L) {
  # 1. simple distance
  count_df_w_dist <- count_df |>
    dplyr::group_by(chr) |>
    dplyr::arrange(chr, coord) |>
    dplyr::mutate(simple_dist = c(NA, diff(coord)))
  curr_group_id <- 0L
  ds <- count_df_w_dist$simple_dist
  group_id <- integer(length = nrow(count_df))
  for (i in seq(1, nrow(count_df))) {
    if (ds[i] > thresh || is.na(ds[i])) {
      curr_group_id <- curr_group_id + 1
    }
    group_id[i] <- curr_group_id
  }
  count_df_w_dist_and_group_string <- count_df_w_dist |>
    dplyr::ungroup() |>
    dplyr::mutate(group_id = group_id) |>
    dplyr::group_by(group_id) |>
    dplyr::mutate(cluster_chr = chr[1],
                  min_cluster_coord = min(coord) - padding,
                  max_cluster_coord = max(coord) + padding) |>
    dplyr::ungroup() |>
    dplyr::mutate(window = paste0(cluster_chr, ":", min_cluster_coord, "-", max_cluster_coord),
                  group_id = NULL, simple_dist = NULL)
  return(count_df_w_dist_and_group_string)
}


#' Construct replicate count table
#'
#' Takes the output of `cluster_loci()` as input; outputs a count matrix for statistical modeling
#'
#' @param clustered_count_df output of `cluster_loci()`, with columns `window`, `umi_count`, and the columns specified by `channel_axes`.
#' @param channel_axes columns whose unique tuples define evidence channels.
#'
#' @returns an integer matrix with evidence channels in the rows and windows in the columns. An entry of the matrix indicates the number of UMIs observed within a given window and channel.
#' @export
#'
#' @examples
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' count_df <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#'  dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & treated & replicate_id %in% 1:2, chr != "chrM") |>
#'  dplyr::select(chr, coord, strand, umi_count, primer_type, replicate_id)
#' clustered_count_df <- cluster_loci(count_df)
#' Y_mat <- construct_replicate_count_table(clustered_count_df)
construct_replicate_count_table <- function(clustered_count_df, channel_axes = c("replicate_id", "primer_type")) {
  x <- lapply(X = channel_axes, FUN = function(col_name) {
    clustered_count_df[[col_name]] |> as.character()
  })
  clustered_count_df$channel_id <- Reduce(f = function(a, b) paste0(a, "_", b), x = x)

  # sum over UMIs within a given (window, channel) pair
  collapsed_count_df <- clustered_count_df |>
    dplyr::group_by(channel_id, window) |>
    dplyr::summarize(umi_count = sum(umi_count), .groups = "drop") |>
    dplyr::arrange(channel_id, window)

  # construct Y_mat
  replicate_idx <- match(x = collapsed_count_df$channel_id, unique(collapsed_count_df$channel_id))
  group_idx <- match(x = collapsed_count_df$window, unique(collapsed_count_df$window))
  Y_mat <- Matrix::sparseMatrix(i = replicate_idx,
                                j = group_idx,
                                x = collapsed_count_df$umi_count) |>
    as.matrix()
  rownames(Y_mat) <- unique(collapsed_count_df$channel_id)
  colnames(Y_mat) <- unique(collapsed_count_df$window)
  return(Y_mat)
}


#' Tune hyperparameters
#'
#' @param Y_mat_trt treated count matrix
#' @param Y_mat_cntrl control count matrix
#' @param c_grid robust hyperparm grid
#' @param lambda_grid lambda grid
#' @param incorporate_occupancy_info a boolean (T/F) indicating whether to incorporate occupancy information into the p-value calculation
#' @param multiplicity_alpha nominal fdr
#' @param max_false_discs maximum false discoveries permitted in the control condition
#' @param annotated_clustered_count_df_trt optional annotated clustered count data frame for the treated condition; if supplied along with `annotated_clustered_count_df_cntrl`, Genovese p-value boosting is used
#' @param annotated_clustered_count_df_cntrl optional annotated clustered count data frame for the control condition; if supplied along with `annotated_clustered_count_df_trt`, Genovese p-value boosting is used
#'
#' @returns a list with elements `selected_params`, `selected_trt_run`, `selected_cntrl_run`, `grid_results`
#' @export
#'
#' @examples
#' # basic
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' count_df_all <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#' dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & replicate_id %in% 1:2, chr != "chrM")
#' Y_mat_trt <- count_df_all |> dplyr::filter(treated) |> cluster_loci() |> construct_replicate_count_table()
#' Y_mat_cntrl <- count_df_all |> dplyr::filter(!treated) |> cluster_loci() |> construct_replicate_count_table()
#' hyperparam_out <- tune_hyperparameters(Y_mat_trt = Y_mat_trt, Y_mat_cntrl = Y_mat_cntrl, c_grid = c(5, 25), lambda_grid = c(5, 50))
#'
#' # with p-value boosting and filtering on homology
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' count_df_all <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#' dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & replicate_id %in% 1:2, chr != "chrM")
#' homology_df <- load_crispritz_output("/Users/timbarry/research_offsite/external/bauer-lab/guideseq_elane/crispritz_CCCCGGCAGAAACGTCCGCG.hg38.targets.txt")
#' n_run_df <- load_n_run_bed("/Users/timbarry/research_offsite/ref_genome_dir/hg38_N_runs_min10.bed")
#' annotated_clustered_count_df_trt <- count_df_all |> dplyr::filter(treated) |> cluster_loci() |>
#'   annotate_clustered_count_df(homology_df = homology_df, n_run_df = n_run_df) |>
#'   dplyr::filter(homology_has_hit)
#' annotated_clustered_count_df_cntrl <- count_df_all |> dplyr::filter(!treated) |> cluster_loci() |>
#'   annotate_clustered_count_df(homology_df = homology_df, n_run_df = n_run_df) |>
#'   dplyr::filter(homology_has_hit)
#' Y_mat_trt <- construct_replicate_count_table(annotated_clustered_count_df_trt)
#' Y_mat_cntrl <- construct_replicate_count_table(annotated_clustered_count_df_cntrl)
#' hyperparam_res <- tune_hyperparameters(Y_mat_trt = Y_mat_trt, Y_mat_cntrl = Y_mat_cntrl,
#'   annotated_clustered_count_df_trt = annotated_clustered_count_df_trt,
#'   annotated_clustered_count_df_cntrl = annotated_clustered_count_df_cntrl)
#'
tune_hyperparameters <- function(Y_mat_trt, Y_mat_cntrl,
                                 c_grid = c(5, 10, 25, 50, 100, 500, 1000),
                                 lambda_grid = c(0, 10, 25, 50, 100),
                                 incorporate_occupancy_info = TRUE,
                                 multiplicity_alpha = 0.5, max_false_discs = 5L,
                                 annotated_clustered_count_df_trt = NULL,
                                 annotated_clustered_count_df_cntrl = NULL,
                                 weight_p_values = TRUE,
                                 lambda_default = 20, prior_strength = "aggressive",
                                 verbose = FALSE) {
  if ((is.null(annotated_clustered_count_df_trt) && !is.null(annotated_clustered_count_df_cntrl)) ||
      (!is.null(annotated_clustered_count_df_trt) && is.null(annotated_clustered_count_df_cntrl))) {
    stop("`annotated_clustered_count_df_trt` and `annotated_clustered_count_df_cntrl` must both be NULL or supplied.")
  }

  # get fitted occupancy models for both conditions
  condition_grid <- c("trt", "cntrl")
  Y_mat_list <- list(trt = Y_mat_trt, cntrl = Y_mat_cntrl)
  annotated_clustered_count_df_list <- list(trt = annotated_clustered_count_df_trt,
                                            cntrl = annotated_clustered_count_df_cntrl)
  occupancy_fit_list <- lapply(X = condition_grid, FUN = function(curr_condition) {
    fit_multirep_guideseq_occupancy(Y_mat = Y_mat_list[[curr_condition]],
                                    incorporate_occupancy_info = incorporate_occupancy_info)
  }) |> setNames(condition_grid)

  # set up grid
  use_occupancy <- sapply(X = occupancy_fit_list, FUN = function(x) {
    x$incorporate_occupancy_info
  })
  if (!all(use_occupancy)) {
    lambda_grid <- lambda_default
    message("Cannot fit occupancy model to both treated and control conditions; fixing lambda to `lambda_default`.")
  }
  grid <- expand.grid(c = c_grid, lambda = lambda_grid, condition = condition_grid)
  fit_grid <- expand.grid(c = c_grid, condition = condition_grid)
  fit_one_count_null <- function(i) {
    print(paste0("Fitting NB model ", i, " of ", nrow(fit_grid)))
    curr_row <- fit_grid[i, , drop = FALSE]
    curr_condition <- as.character(curr_row$condition)
    curr_c <- curr_row$c[[1]]
    mu_theta_hat <- fit_multirep_guideseq_count_null(Y_mat = Y_mat_list[[curr_condition]],
                                                     c_tukey_beta = curr_c,
                                                     c_tukey_sigma = curr_c,
                                                     robust_fit = TRUE)
    list(params = curr_row, mu_theta_hat = mu_theta_hat)
  }

  count_fit_list <- lapply(X = seq_len(nrow(fit_grid)), FUN = fit_one_count_null)
  count_fit_names <- sapply(count_fit_list, FUN = function(curr_fit) {
    paste(as.character(curr_fit$params$condition), curr_fit$params$c[[1]], sep = "_")
  })
  names(count_fit_list) <- count_fit_names

  score_one_grid_row <- function(i) {
    print(paste0("Scoring ", i, " of ", nrow(grid)))
    curr_row <- grid[i, , drop = FALSE]
    curr_condition <- as.character(curr_row$condition)
    curr_c <- curr_row$c[[1]]
    curr_lambda <- curr_row$lambda[[1]]
    mu_theta_hat <- count_fit_list[[paste(curr_condition, curr_c, sep = "_")]]$mu_theta_hat
    if (curr_condition == "trt") {
      curr_Y_mat <- Y_mat_trt
    } else {
      curr_Y_mat <- Y_mat_cntrl
    }
    annotated_clustered_count_df <- annotated_clustered_count_df_list[[curr_condition]]
    fit_res <- score_multirep_guideseq_fit(Y_mat = curr_Y_mat,
                                           occupancy_fit = occupancy_fit_list[[curr_condition]],
                                           mu_theta_hat_mat = mu_theta_hat,
                                           lambda = curr_lambda,
                                           multiplicity_alpha = multiplicity_alpha,
                                           annotated_clustered_count_df = annotated_clustered_count_df)
    if (weight_p_values && !is.null(annotated_clustered_count_df_trt) && !is.null(annotated_clustered_count_df_cntrl)) {
      fit_res$res_df <- fit_res$res_df |> boost_p_values_genovese_cfd(multiplicity_alpha = multiplicity_alpha, prior_strength = prior_strength)
    }
    list(params = curr_row, res = fit_res)
  }

  grid_results <- lapply(X = seq_len(nrow(grid)), FUN = score_one_grid_row)
  summary_df <- lapply(X = grid_results, FUN = function(curr_res) {
    curr_res$params |>
      dplyr::mutate(n_discoveries = sum(curr_res$res$res_df$nominated_window))
  }) |> data.table::rbindlist() |>
    tidyr::pivot_wider(names_from = "condition", values_from = n_discoveries)
  if (any(summary_df$cntrl <= max_false_discs)) {
    selected_params <- summary_df |>
      dplyr::filter(cntrl <= max_false_discs) |>
      dplyr::arrange(dplyr::desc(trt), dplyr::desc(c), lambda) |>
      dplyr::slice(1)
    trt_idx <- sapply(grid_results, FUN = function(curr_res) {
      curr_res$params$c == selected_params$c && curr_res$params$lambda == selected_params$lambda && curr_res$params$condition == "trt"
    }) |> which()
    cntrl_idx <- sapply(grid_results, FUN = function(curr_res) {
      curr_res$params$c == selected_params$c && curr_res$params$lambda == selected_params$lambda && curr_res$params$condition == "cntrl"
    }) |> which()
    selected_trt_run <- grid_results[[trt_idx]]$res
    selected_cntrl_run <- grid_results[[cntrl_idx]]$res
  } else {
    selected_params <- NA
    selected_trt_run <- NA
    selected_cntrl_run <- NA
  }
  ret <- list(selected_params = selected_params,
              selected_trt_run = selected_trt_run,
              selected_cntrl_run = selected_cntrl_run,
              grid_results = grid_results, summary_df = summary_df)
  return(ret)
}


#' Load CRISPRitz output
#'
#' @param targets_file_path file path to the targets.txt file outputted by CRISPRitz
#'
#' @returns a data frame containing the CRISPRitz target output
#' @export
#'
#' @examples
#' homology_df <- load_crispritz_output("/Users/timbarry/research_offsite/external/bauer-lab/guideseq_elane/crispritz_CCCCGGCAGAAACGTCCGCG.hg38.targets.txt")
load_crispritz_output <- function(targets_file_path) {
  df <- readr::read_delim(file = targets_file_path) |>
    dplyr::rename("bulge_type" = "#Bulge type", "n_mismatches" = "Mismatches",
                  "n_bulges" = "Bulge Size", "n_total_changes" = "Total",
                  "gRNA" = "crRNA", "dna" = "DNA", "chromosome" = "Chromosome",
                  "posit" = "Position", "cluster_posit" = "Cluster Position",
                  "strand" = "Direction")
  protospacer_width <- nchar(gsub(pattern = "-", replacement = "", x = df$dna)) - 3L
  df$protospacer_width <- protospacer_width
  return(df)
}


#' Load N run bed file
#'
#' Load bed file storing the positions of N-runs in the reference genome.
#'
#' @param n_run_bed_file_path path to N-run bed file
#'
#' @returns a data frame storing the positions of the N-runs
#' @export
#'
#' @examples
#' n_run_bed_file_path <- "/Users/timbarry/research_offsite/ref_genome_dir/hg38_N_runs_min10.bed"
#' n_run_df <- load_n_run_bed(n_run_bed_file_path)
load_n_run_bed <- function(n_run_bed_file_path) {
  df <- readr::read_delim(file = n_run_bed_file_path, col_names = FALSE)
  colnames(df) <- c("chromosome", "start", "end", "feature", "score", "strand")
  df |> dplyr::mutate(start = start + 1L)
}


#' Load ENCODE blacklist bed
#'
#' @param encode_blacklist_file_path filepath to ENCODE blacklist
#'
#' @returns a data frame storing the positions of the blacklisted regions
#' @export
#'
#' @examples
#' encode_blacklist_file_path <- "/Users/timbarry/research_offsite/ref_genome_dir/hg38-blacklist.v2.bed"
load_encode_blacklist_bed <- function(encode_blacklist_file_path) {
  df <- readr::read_delim(file = encode_blacklist_file_path, col_names = FALSE)
  colnames(df) <- c("chromosome", "start", "end", "classification")
  df |> dplyr::mutate(start = start + 1L)
}

#' Annotate clustered count df with homology
#'
#' Add columns:
#' - `homology_has_hit` (indicating whether there is a homology hit inside the window)
#' - `overlaps_n_run` (indicating whether the window overlaps an N-run in the reference genome)
#' - `overlaps_encode_blacklist` (indicating whether the window overlaps an ENCODE blacklist region)
#' - `homology_n_mismatches` (indicating the number of mismatches between aligned spacer and protospacer sequence)
#' - `homology_n_bulges` (indicating number of bulges between aligned spacer and protospacer sequence)
#' - `homology_bulge_type` (indicating the CRISPRitz bulge type)
#' - `homology_cfd` (indicating the CFD-like homology score)
#' - `homology_chromosome` (chromosome of the aligned protospacer)
#' - `homology_strand` (indicating whether the protospacer is on the plus or minus strand)
#' - `homology_dna` (aligned protospacer sequence)
#' - `homology_gRNA` (aligned spacer sequence)
#' - `homology_cut_start` and `homology_cut_end` (predicted cut-site coordinates)
#' - `homology_modal_base_cut_distance` (distance between the modal base and predicted cut site)
#' - `homology_alignment_score` (combined homology and cut-site-distance score)
#'
#' @param clustered_count_df output of `cluster_loci()`
#' @param homology_df optional output of `load_crispritz_output()`; if supplied, windows are annotated for overlap with CRISPRitz hits
#' @param n_run_df optional output of `load_n_run_bed()`; if supplied, windows are annotated for overlap with N-runs
#' @param encode_blacklist_df optional output of `load_encode_blacklist_bed()`; if supplied, windows are annotated for overlap with ENCODE blacklist regions
#'
#' @examples
#' homology_df <- load_crispritz_output("/Users/timbarry/research_offsite/external/bauer-lab/guideseq_bcl11a/1620_crispritz_spRY_bcl11a_windows.hg38.targets.txt")
#' n_run_df <- load_n_run_bed("/Users/timbarry/research_offsite/ref_genome_dir/hg38_N_runs_min10.bed")
#' encode_blacklist_df <- load_encode_blacklist_bed("/Users/timbarry/research_offsite/ref_genome_dir/hg38-blacklist.v2.bed")
#' bcl11a_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_bcl11a/")
#' clustered_count_df <- readRDS(paste0(bcl11a_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#' dplyr::filter(grna == "1620" & treated & chr != "chrM") |>
#' dplyr::select(chr, coord, strand, umi_count, primer_type, replicate_id) |>
#' cluster_loci()
#' annotated_clustered_count_df <- annotate_clustered_count_df(clustered_count_df = clustered_count_df,
#'   homology_df = homology_df, n_run_df = n_run_df, encode_blacklist_df = encode_blacklist_df)
#' @export
annotate_clustered_count_df <- function(clustered_count_df, homology_df = NULL, n_run_df = NULL, encode_blacklist_df = NULL) {
  # compute window df by computing a summary over clustered_count_df
  window_df <- clustered_count_df |>
    dplyr::group_by(window, coord) |>
    dplyr::summarize(umi_count_per_base = sum(umi_count),
                     chr = chr[1],
                     max_cluster_coord = max_cluster_coord[1],
                     min_cluster_coord = min_cluster_coord[1]) |>
    dplyr::summarize(modal_base = coord[umi_count_per_base == max(umi_count_per_base)][1],
                     window_n_occupied_bases = length(unique(coord)),
                     window_width = max_cluster_coord[1] - min_cluster_coord[1] + 1L,
                     chr = chr[1],
                     min_cluster_coord = min_cluster_coord[1],
                     max_cluster_coord = max_cluster_coord[1])

  # construct window gr
  window_gr <- GenomicRanges::GRanges(
    seqnames = window_df$chr,
    ranges = IRanges::IRanges(start = window_df$min_cluster_coord, end = window_df$max_cluster_coord)
  )

  append_simple_overlap <- function(window_df, window_gr, df_in, colname) {
    gr_in <- GenomicRanges::GRanges(
      seqnames = df_in$chromosome,
      ranges = IRanges::IRanges(start = df_in$start, end = df_in$end)
    )
    seqlevels <- union(GenomeInfoDb::seqlevels(window_gr), GenomeInfoDb::seqlevels(gr_in))
    GenomeInfoDb::seqlevels(window_gr) <- seqlevels
    GenomeInfoDb::seqlevels(gr_in) <- seqlevels
    hits <- GenomicRanges::findOverlaps(query = window_gr, subject = gr_in,
                                        ignore.strand = TRUE)
    window_df[[colname]] <- FALSE
    window_df[[colname]][S4Vectors::queryHits(hits)] <- TRUE
    return(window_df)
  }
  if (!is.null(n_run_df)) {
    window_df <- append_simple_overlap(window_df = window_df,
                                       window_gr = window_gr,
                                       df_in = n_run_df,
                                       colname = "overlaps_n_run")
  }
  if (!is.null(encode_blacklist_df)) {
    window_df <- append_simple_overlap(window_df = window_df,
                                       window_gr = window_gr,
                                       df_in = encode_blacklist_df,
                                       colname = "overlaps_encode_blacklist")
  }

  # add homology information if homology_df provided
  if (!is.null(homology_df)) {
    # construct homology gr object
    homology_gr <- GenomicRanges::GRanges(
      seqnames = homology_df$chromosome,
      ranges = IRanges::IRanges(
        start = homology_df$posit + 1L + ifelse(homology_df$strand == "-", 3L, 0L),
        width = homology_df$protospacer_width),
      strand = homology_df$strand)

    # find overlaps between homology gr and window gr; retain alignments overlapping a GUIDE-seq window; append cfd score
    crispritz_hits <- GenomicRanges::findOverlaps(query = window_gr,
                                                  subject = homology_gr,
                                                  ignore.strand = TRUE)
    if (length(crispritz_hits) == 0L) {
      window_df <- window_df |>
        dplyr::mutate(homology_has_hit = FALSE,
                      homology_bulge_type = NA_character_, homology_gRNA = NA_character_,
                      homology_dna = NA_character_, homology_chromosome = NA_character_,
                      homology_strand = NA_character_, homology_n_mismatches = NA_real_,
                      homology_n_bulges = NA_real_, homology_cut_start = NA_real_,
                      homology_cut_end = NA_real_, homology_cfd = NA_real_,
                      homology_modal_base_cut_distance = NA_real_, homology_alignment_score = NA_real_)
    } else {
      crispritz_subject_hits <- S4Vectors::subjectHits(crispritz_hits)
      crispritz_query_hits <- S4Vectors::queryHits(crispritz_hits)
      homology_df_sub <- homology_df[crispritz_subject_hits,] |>
        dplyr::mutate(homology_cut_start = posit + 1L + ifelse(strand == "+", protospacer_width - 4L, 5L),
                      homology_cut_end = posit + 1L + ifelse(strand == "+", protospacer_width - 3L, 6L)) |>
        dplyr::mutate(homology_cfd = calculate_cfd_score(homology_dna = dna, homology_gRNA = gRNA, homology_has_hit = TRUE))

      # for each window (containing a CRISPRitz hit), score each alignment and return the best
      window_idxs_with_hit <- unique(crispritz_query_hits)
      homology_df <- lapply(X = window_idxs_with_hit, FUN = function(i) {
        curr_window <- window_df[i,]
        curr_crispritz_candidates <- homology_df_sub[crispritz_query_hits == i,] |>
          dplyr::mutate(homology_modal_base_cut_distance = pmin(abs(curr_window$modal_base - homology_cut_start),
                                                                abs(curr_window$modal_base - homology_cut_end))) |>
          dplyr::mutate(homology_alignment_score = compute_alignment_scores(cfds = homology_cfd, distances = homology_modal_base_cut_distance))
        # find the alignment with the best score, tie-breaking by distance
        best_alignment <- curr_crispritz_candidates |>
          dplyr::arrange(dplyr::desc(homology_alignment_score), homology_modal_base_cut_distance, dplyr::desc(homology_cfd)) |>
          dplyr::slice(1L) |>
          dplyr::select(bulge_type, gRNA, dna, chromosome, strand, n_mismatches, n_bulges,
                        homology_cut_start, homology_cut_end, homology_cfd, homology_modal_base_cut_distance,
                        homology_alignment_score) |>
          dplyr::rename(homology_bulge_type = bulge_type, homology_gRNA = gRNA, homology_dna = dna,
                        homology_chromosome = chromosome, homology_strand = strand, homology_n_mismatches = n_mismatches,
                        homology_n_bulges = n_bulges) |>
          dplyr::mutate(window = curr_window$window)
      }) |> data.table::rbindlist()
      window_df <- dplyr::left_join(window_df, homology_df, by = "window")
    }
  }

  # prepare final output
  window_df <- window_df |> dplyr::select(-min_cluster_coord, -max_cluster_coord, -chr)
  clustered_count_df <- dplyr::left_join(x = clustered_count_df, y = window_df, by = "window")
  if (!is.null(homology_df)) {
    clustered_count_df <- clustered_count_df |>
      dplyr::mutate(homology_has_hit = ifelse(is.na(homology_bulge_type), FALSE, TRUE))
  }

  return(clustered_count_df)
}
