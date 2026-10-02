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
#' @param tau baseline normalized weight for zero homology scores
#' @param gamma exponential distance-decay coefficient; annotations should be generated with the same value
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
#' # future::plan(future::multisession, workers = 4)
#' future::plan(future::sequential)
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
                                 lambda_default = 20, tau = 0.1, gamma = log(20)/7,
                                 verbose = FALSE) {
  if ((is.null(annotated_clustered_count_df_trt) && !is.null(annotated_clustered_count_df_cntrl)) ||
      (!is.null(annotated_clustered_count_df_trt) && is.null(annotated_clustered_count_df_cntrl))) {
    stop("`annotated_clustered_count_df_trt` and `annotated_clustered_count_df_cntrl` must both be NULL or supplied.")
  }

  ###########################################################
  # PART 1: FIT OCCUPANCY AND COUNT MODELS TO BOTH CONDITIONS
  ###########################################################
  # fit occupancy models to both conditions
  condition_grid <- c("trt", "cntrl")
  Y_mat_list <- list(trt = Y_mat_trt, cntrl = Y_mat_cntrl)
  annotated_clustered_count_df_list <- list(trt = annotated_clustered_count_df_trt,
                                            cntrl = annotated_clustered_count_df_cntrl)
  occupancy_fit_list <- lapply(X = condition_grid, FUN = function(curr_condition) {
    fit_multirep_guideseq_occupancy(Y_mat = Y_mat_list[[curr_condition]],
                                    incorporate_occupancy_info = incorporate_occupancy_info)
  }) |> setNames(condition_grid)

  # determine whether to use occupancy
  use_occupancy <- sapply(X = occupancy_fit_list, FUN = function(x) {
    x$incorporate_occupancy_info
  })
  if (!all(use_occupancy)) {
    lambda_grid <- lambda_default
    message("Cannot fit occupancy model to both treated and control conditions; fixing lambda to `lambda_default`.")
  }

  # fit the NB models, iterating over c_grid
  nb_model_fits <- future.apply::future_lapply(X = c_grid, FUN = function(curr_c) {
    message("Running NB fits for c = ", curr_c)
    lapply(X = condition_grid, FUN = function(curr_condition) {
      mu_theta_hat <- fit_multirep_guideseq_count_null(Y_mat = Y_mat_list[[curr_condition]],
                                                       c_tukey_beta = curr_c,
                                                       c_tukey_sigma = curr_c,
                                                       robust_fit = TRUE)
      return(mu_theta_hat)
    }) |> setNames(condition_grid)
  }) |> setNames(c_grid)

  ####################################################################
  # PART 2: COMPUTE TEST STATISTICS (WHICH DO NOT DEPEND NB MODEL FIT)
  ####################################################################
  test_stat_list <- lapply(X = condition_grid, FUN = function(condition) {
    # get the count matrix and occupancy fit
    Y_mat <- Y_mat_list[[condition]]
    occupancy_fit <- occupancy_fit_list[[condition]]
    total_umi_counts <- colSums(Y_mat)
    occupancy_counts <- colSums(occupancy_fit$X)
    pattern_log_pi_sum <- test_stats <- test_stats_by_lambda <- NULL
    # (i) for incorporate occupancy info or not, compute test stats, then compute max_needed
    if (!occupancy_fit$incorporate_occupancy_info) { # occupancy-blind model
      test_stats <- total_umi_counts - occupancy_counts
      max_needed <- max(test_stats)
    } else { # occupancy-aware model
      log_pi_hat <- log(occupancy_fit$pi_hat)
      window_log_pi_sum <- as.numeric(crossprod(log_pi_hat, occupancy_fit$X))
      pattern_log_pi_sum <- as.numeric(occupancy_fit$Omega %*% log_pi_hat)
      # compute the test statistics over lambda
      test_stats_by_lambda <- lapply(X = lambda_grid, FUN = function(lambda) {
        test_stats <- (total_umi_counts - occupancy_counts) - lambda * window_log_pi_sum
      }) |> setNames(lambda_grid)
      max_test_stat_by_lambda <- sapply(X = test_stats_by_lambda, FUN = max)
      max_needed <- max(0L, ceiling(max_test_stat_by_lambda + lambda_grid * max(pattern_log_pi_sum)))
    }
    return(list(test_stats_by_lambda = test_stats_by_lambda,
                test_stats = test_stats, max_needed = max_needed,
                pattern_log_pi_sum = pattern_log_pi_sum))
  }) |> setNames(condition_grid)

  ##########################
  # PART 3: COMPUTE P-VALUES
  ##########################
  score_model_for_given_c <- function(c) {
    message(paste0("Scoring c = ", c, ", lambda = ", paste0(lambda_grid,collapse = ", ")))
    # iterate over conditions
    p_vals_by_condition <- lapply(X = condition_grid, FUN = function(condition) {
      mu_theta_hat_mat <- nb_model_fits[[as.character(c)]][[condition]]
      max_needed <- test_stat_list[[condition]]$max_needed
      occupancy_fit <- occupancy_fit_list[[condition]]
      # get the list of partial convolutions
      right_tail_prob_list <- get_right_tail_prob_list(mu_theta_hat_mat = mu_theta_hat_mat,
                                                       max_needed = max_needed,
                                                       Omega = occupancy_fit$Omega)
      p_vals <- p_vals_per_lambda <- NULL
      if (!occupancy_fit$incorporate_occupancy_info) { # not incorporating occupancy info
        occupancy_pattern_map <- occupancy_fit$occupancy_pattern_map
        test_stats <- test_stat_list[[condition]]$test_stats
        p_vals <- numeric(length(occupancy_pattern_map))
        for (i in seq_along(right_tail_prob_list)) {
          idxs <- which(occupancy_pattern_map == i)
          p_vals[idxs] <- get_p_values_given_test_stats_prob_vector(
            test_stat_v_in = test_stats[idxs],
            right_tail_prob_v_in = right_tail_prob_list[[i]]
          )
          p_vals <- pmin(1, p_vals)
        }
      } else {
        # loop over lambda
        p_vals_per_lambda <- lapply(X = lambda_grid, FUN = function(lambda) {
          test_stats <- test_stat_list[[condition]]$test_stats_by_lambda[[as.character(lambda)]]
          pattern_log_pi_sum <- test_stat_list[[condition]]$pattern_log_pi_sum
          l <- sapply(X = seq_len(nrow(occupancy_fit$Omega)), FUN = function(i) {
            sum_start <- ceiling(test_stats + lambda * pattern_log_pi_sum[i])
            nb_piece <- get_p_values_given_test_stats_prob_vector(
              test_stat_v_in = sum_start,
              right_tail_prob_v_in = right_tail_prob_list[[i]]
            )
            occupancy_fit$tbp_pattern_df$pmf[i] * nb_piece
          }, simplify = FALSE)
          p_vals <- pmin(1, Reduce(f = "+", x = l))
        }) |> setNames(lambda_grid)
      }
      return(list(p_vals_per_lambda = p_vals_per_lambda, p_vals = p_vals))
    }) |> setNames(condition_grid)
  }
  p_vals_by_c <- future.apply::future_lapply(X = c_grid, FUN = score_model_for_given_c) |>
    setNames(c_grid)

  ########################
  # PART 4: PREPARE OUTPUT
  ########################
  message("Collating results and preparing outout.")
  grid <- expand.grid(c = c_grid, lambda = lambda_grid, condition = condition_grid)
  grid_results <- lapply(X = seq_len(nrow(grid)), FUN = function(i) {
    curr_row <- grid[i,,drop = FALSE]
    curr_condition <- as.character(curr_row$condition)
    curr_c <- as.character(curr_row$c)
    curr_lambda <- as.character(curr_row$lambda)
    Y_mat <- Y_mat_list[[curr_condition]]
    occupancy_fit <- occupancy_fit_list[[curr_condition]]
    curr_p_vals <- p_vals_by_c[[curr_c]][[curr_condition]]
    curr_test_stats <- test_stat_list[[curr_condition]]
    if (occupancy_fit$incorporate_occupancy_info) {
      p_vals <- curr_p_vals$p_vals_per_lambda[[curr_lambda]]
      test_stats <- curr_test_stats$test_stats_by_lambda[[curr_lambda]]
    } else {
      p_vals <- curr_p_vals$p_vals
      test_stats <- curr_test_stats$test_stats
    }

    res_df <- data.frame(window = colnames(Y_mat),
                         p_value = p_vals,
                         test_stat = test_stats,
                         nominated_window = p.adjust(p = p_vals, method = "BH") < multiplicity_alpha,
                         umi_count = colSums(Y_mat),
                         lambda = if (occupancy_fit$incorporate_occupancy_info) curr_row$lambda else NA,
                         occupancy_pattern = occupancy_fit$col_keys) |> dplyr::arrange(p_value)
    rownames(res_df) <- NULL
    annotated_clustered_count_df <- annotated_clustered_count_df_list[[curr_condition]]
    if (!is.null(annotated_clustered_count_df)) {
      right_df <- annotated_clustered_count_df |>
        dplyr::select(window, dplyr::starts_with(c("homology", "window", "overlaps"))) |>
        dplyr::filter(window %in% res_df$window) |>
        dplyr::distinct()
      res_df <- dplyr::left_join(res_df, right_df, by = "window")
      if (weight_p_values) {
        res_df <- boost_p_values_genovese_cfd(augmented_result_df = res_df,
                                             multiplicity_alpha = multiplicity_alpha,
                                             tau = tau, gamma = gamma)
      }
    }
    ests_list <- list(mu_theta_hat_mat = nb_model_fits[[curr_c]][[curr_condition]])
    if (occupancy_fit$incorporate_occupancy_info) ests_list$pi_hat <- occupancy_fit$pi_hat
    list(params = curr_row, res = list(res_df = res_df, ests_list = ests_list))
  })

  summary_df <- lapply(X = grid_results, FUN = function(curr_res) {
    curr_res$params |>
      dplyr::mutate(n_discoveries = sum(curr_res$res$res_df$nominated_window))
  }) |> data.table::rbindlist() |>
    tidyr::pivot_wider(names_from = "condition", values_from = n_discoveries)
  if (any(summary_df$cntrl <= max_false_discs)) {
    selected_params <- summary_df |>
      dplyr::filter(cntrl <= max_false_discs) |>
      dplyr::arrange(cntrl, dplyr::desc(trt), dplyr::desc(c), lambda) |>
      dplyr::arrange(dplyr::desc(trt), cntrl, dplyr::desc(c), lambda) |>
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


#' Score multirep GUIDE-seq fit
#'
#' @param Y_mat Y matrix
#' @param occupancy_fit fitted occupancy model
#' @param mu_theta_hat_mat matrix of fitted NB model parameters
#' @param lambda_grid grid of lambda values over which to compute p-value
#'
#' @returns
#' @export
#'
#' @examples
score_multirep_guideseq_fit <- function(Y_mat, occupancy_fit, mu_theta_hat_mat, lambda_grid) {
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
                                                     max_needed = max(test_stats), Omega = Omega)
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

