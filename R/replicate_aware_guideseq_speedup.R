#' Tune hyperparameters
#'
#' @param Y_mat_trt treated count matrix
#' @param Y_mat_cntrl control count matrix
#' @param c_grid robust hyperparm grid
#' @param lambda_grid lambda grid
#' @param incorporate_occupancy_info a boolean (T/F) indicating whether to incorporate occupancy information into the p-value calculation
#' @param multiplicity_alpha nominal fdr
#' @param max_false_discs maximum false discoveries permitted in the control condition
#' @param annotated_clustered_count_df_trt annotated clustered count data frame for the treated condition
#' @param annotated_clustered_count_df_cntrl annotated clustered count data frame for the control condition
#' @param tau baseline normalized weight for zero homology scores
#' @param verbose whether to print progress messages
#'
#' @returns A list containing `result_list` (skinny window-level tables indexed by
#'   condition, c, and lambda), `nb_model_fits` (parameter matrices indexed by
#'   condition and c), `summary_df`, `selected_params`, `selected_trt_run`, and
#'   `selected_cntrl_run`. Selected runs contain full annotations sorted by
#'   p-value and their fitted parameters. If no combination meets the
#'   control-discovery limit, the selected outputs are NA.
#' @export
#'
#' @examples
#' # with p-value boosting and filtering on homology
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' count_df_all <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#' dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & replicate_id %in% 1:2, chr != "chrM")
#' homology_df <- load_crispritz_output("/Users/timbarry/research_offsite/external/bauer-lab/guideseq_elane/crispritz_CCCCGGCAGAAACGTCCGCG.hg38.targets.txt")
#' n_run_df <- load_n_run_bed("/Users/timbarry/research_offsite/ref_genome_dir/hg38_N_runs_min10.bed")
#' annotated_clustered_count_df_trt <- count_df_all |> dplyr::filter(treated) |> cluster_loci() |>
#'   annotate_clustered_count_df(homology_df = homology_df, n_run_df = n_run_df)
#' annotated_clustered_count_df_cntrl <- count_df_all |> dplyr::filter(!treated) |> cluster_loci() |>
#'   annotate_clustered_count_df(homology_df = homology_df, n_run_df = n_run_df)
#' Y_mat_trt <- construct_replicate_count_table(annotated_clustered_count_df_trt)
#' Y_mat_cntrl <- construct_replicate_count_table(annotated_clustered_count_df_cntrl)
#' hyperparam_res <- tune_hyperparameters(Y_mat_trt = Y_mat_trt, Y_mat_cntrl = Y_mat_cntrl,
#'   annotated_clustered_count_df_trt = annotated_clustered_count_df_trt,
#'   annotated_clustered_count_df_cntrl = annotated_clustered_count_df_cntrl)
tune_hyperparameters <- function(Y_mat_trt, Y_mat_cntrl,
                                 annotated_clustered_count_df_trt,
                                 annotated_clustered_count_df_cntrl,
                                 c_grid = c(5, 10, 25, 50, 100, 500, 1000),
                                 lambda_grid = c(0, 10, 25, 50, 100),
                                 incorporate_occupancy_info = TRUE,
                                 multiplicity_alpha = 0.5, max_false_discs = 5L,
                                 weight_p_values = TRUE,
                                 lambda_default = 20, tau = 0.1,
                                 verbose = FALSE) {
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

  # fit the robust NB models, iterating over c_grid
  nb_model_fits <- lapply(X = condition_grid, FUN = function(curr_condition) {
    if (verbose) message("Fitting NB models for ", curr_condition)
    fit_multirep_guideseq_count_null(Y_mat = Y_mat_list[[curr_condition]], c_grid = c_grid)
  }) |> setNames(condition_grid)

  ##############################################################
  # PART 2: GROUP OBSERVATIONS AND COMPUTE TEST STATISTICS
  ##############################################################
  observation_group_list <- lapply(X = condition_grid, FUN = function(condition) {
    umi_counts <- colSums(Y_mat_list[[condition]])
    occupancy_patterns <- occupancy_fit_list[[condition]]$col_keys
    occupancy_pattern_map <- occupancy_fit_list[[condition]]$occupancy_pattern_map
    windows <- names(umi_counts)
    names(occupancy_pattern_map) <- names(umi_counts) <- names(occupancy_patterns) <- NULL
    window_to_observation_mapping_df <- data.frame(window = windows,
                                                   umi_count = umi_counts,
                                                   occupancy_pattern = occupancy_patterns,
                                                   occupancy_pattern_map = occupancy_pattern_map)
    unique_observation_df <- window_to_observation_mapping_df |> dplyr::select(-window) |>
      dplyr::distinct(umi_count, occupancy_pattern, occupancy_pattern_map)
    out <- list(window_to_observation_mapping_df = window_to_observation_mapping_df,
                unique_observation_df = unique_observation_df)
    return(out)
  }) |> setNames(condition_grid)

  # compute test statistics for the unique observations
  test_stat_list <- lapply(X = condition_grid, FUN = function(condition) {
    unique_observation_df <- observation_group_list[[condition]]$unique_observation_df
    occupancy_fit <- occupancy_fit_list[[condition]]
    X_unique <- strsplit(x = unique_observation_df$occupancy_pattern, split = "") |>
      sapply(as.integer)
    occupancy_count <- colSums(X_unique)
    pattern_log_pi_sum <- test_stats <- test_stats_by_lambda <- NULL
    # (i) for incorporate occupancy info or not, compute test stats, then compute max_needed
    if (!occupancy_fit$incorporate_occupancy_info) { # occupancy-blind model
      test_stats <- unique_observation_df$umi_count - occupancy_count
      max_needed <- max(test_stats)
    } else { # occupancy-aware model
      log_pi_hat <- log(occupancy_fit$pi_hat)
      group_log_pi_sum <- as.numeric(crossprod(log_pi_hat, X_unique))
      pattern_log_pi_sum <- as.numeric(occupancy_fit$Omega %*% log_pi_hat)
      # compute the test statistics over lambda
      test_stats_by_lambda <- lapply(X = lambda_grid, FUN = function(lambda) {
        test_stats <- (unique_observation_df$umi_count - occupancy_count) - lambda * group_log_pi_sum
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
    if (verbose) message("Scoring c = ", c)
    # iterate over conditions
    p_vals_by_condition <- lapply(X = condition_grid, FUN = function(condition) {
      mu_theta_hat_mat <- nb_model_fits[[condition]][[as.character(c)]]
      max_needed <- test_stat_list[[condition]]$max_needed
      occupancy_fit <- occupancy_fit_list[[condition]]
      # get the list of partial convolutions
      right_tail_prob_list <- get_right_tail_prob_list(mu_theta_hat_mat = mu_theta_hat_mat,
                                                       max_needed = max_needed,
                                                       Omega = occupancy_fit$Omega)
      p_vals <- p_vals_per_lambda <- NULL
      if (!occupancy_fit$incorporate_occupancy_info) { # not incorporating occupancy info
        test_stats <- test_stat_list[[condition]]$test_stats
        occupancy_pattern_map <- observation_group_list[[condition]]$unique_observation_df$occupancy_pattern_map
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
  p_vals_by_c <- lapply(X = c_grid, FUN = score_model_for_given_c) |> setNames(c_grid)

  ##################################
  # PART 4: PREPARE WINDOW RESULTS
  ##################################
  if (verbose) message("Preparing window results")
  # construct the starting result df
  result_dfs <- lapply(X = condition_grid, FUN = function(condition) {
    # starting point: (window, umi_count, occupancy_pattern) df
    result_df <- observation_group_list[[condition]]$window_to_observation_mapping_df |>
      dplyr::select(-occupancy_pattern_map)
    # append the group_id of each window
    unique_group_df <- observation_group_list[[condition]]$unique_observation_df |>
      dplyr::select(-occupancy_pattern_map)
    unique_group_df$group_id <- seq_len(nrow(unique_group_df))
    result_df <- dplyr::left_join(x = result_df, y = unique_group_df,
                                  by = c("umi_count", "occupancy_pattern"))
    # join with annotated df
    annotated_clustered_count_df <- annotated_clustered_count_df_list[[condition]]
    annotation_df <- annotated_clustered_count_df |>
      dplyr::select(window, dplyr::starts_with(c("homology", "window", "overlaps"))) |>
      dplyr::distinct()
    result_df <- dplyr::left_join(x = result_df, y = annotation_df, by = "window")
  }) |> setNames(condition_grid)

  # expand p-values to windows and apply BH and homology weighting
  result_list <- lapply(X = condition_grid, FUN = function(condition) {
    result_df_skinny <- result_dfs[[condition]] |> dplyr::select(window, group_id, homology_alignment_score)
    incorporate_occupancy_info <- occupancy_fit_list[[condition]]$incorporate_occupancy_info
    results_over_c_lambda <- lapply(X = c_grid, FUN = function(c) {
      lapply(X = lambda_grid, FUN = function(lambda) {
        # extract the relevant p-values for this (c, lambda)
        p_vals_per_lambda <- p_vals_by_c[[as.character(c)]][[condition]]
        if (incorporate_occupancy_info) {
          p_vals <- p_vals_per_lambda$p_vals_per_lambda[[as.character(lambda)]]
        } else {
          p_vals <- p_vals_per_lambda$p_vals
        }
        result_df_skinny <- result_df_skinny |> dplyr::mutate(p_value = p_vals[group_id])
        if (weight_p_values) {
          result_df_skinny <- boost_p_values_genovese_cfd(augmented_result_df = result_df_skinny,
                                                          multiplicity_alpha = multiplicity_alpha, tau = tau)
        } else {
          result_df_skinny$nominated_window <- p.adjust(p = result_df_skinny$p_value, method = "BH") < multiplicity_alpha
        }
        result_df_skinny
      }) |> setNames(lambda_grid)
    }) |> setNames(c_grid)
  }) |> setNames(condition_grid)


  #################################
  # PART 5: SUMMARIZE AND SELECT
  #################################
  summary_df <- expand.grid(c = c_grid, lambda = lambda_grid, KEEP.OUT.ATTRS = FALSE) |> dplyr::as_tibble()
  for (condition in condition_grid) {
    summary_df[[condition]] <- vapply(seq_len(nrow(summary_df)), function(i) {
      curr_c <- as.character(summary_df$c[i])
      curr_lambda <- as.character(summary_df$lambda[i])
      sum(result_list[[condition]][[curr_c]][[curr_lambda]]$nominated_window)
    }, integer(1))
  }
  if (any(summary_df$cntrl <= max_false_discs)) {
    selected_params <- summary_df |>
      dplyr::filter(cntrl <= max_false_discs) |>
      dplyr::arrange(dplyr::desc(trt), cntrl, dplyr::desc(c), lambda) |>
      dplyr::slice(1)

    # join metadata only for the selected runs
    selected_c <- as.character(selected_params$c)
    selected_lambda <- as.character(selected_params$lambda)
    selected_runs <- lapply(X = condition_grid, FUN = function(condition) {
      res_df <- result_list[[condition]][[selected_c]][[selected_lambda]]
      occupancy_fit <- occupancy_fit_list[[condition]]
      ests_list <- list(mu_theta_hat_mat = nb_model_fits[[condition]][[selected_c]])
      if (occupancy_fit$incorporate_occupancy_info) ests_list$pi_hat <- occupancy_fit$pi_hat
      res_df$lambda <- if (occupancy_fit$incorporate_occupancy_info) selected_params$lambda else NA
      metadata_df <- result_dfs[[condition]] |>
        dplyr::select(-group_id, -homology_alignment_score)
      res_df <- res_df |>
        dplyr::select(-group_id) |>
        dplyr::left_join(metadata_df, by = "window") |>
        dplyr::relocate(window, p_value, nominated_window, umi_count, lambda, occupancy_pattern) |>
        dplyr::arrange(p_value)
      list(res_df = res_df, ests_list = ests_list)
    }) |> setNames(condition_grid)
    selected_trt_run <- selected_runs$trt
    selected_cntrl_run <- selected_runs$cntrl
  } else {
    selected_params <- NA
    selected_trt_run <- NA
    selected_cntrl_run <- NA
  }
  ret <- list(selected_params = selected_params,
              selected_trt_run = selected_trt_run,
              selected_cntrl_run = selected_cntrl_run,
              result_list = result_list, nb_model_fits = nb_model_fits,
              summary_df = summary_df)
  return(ret)
}


get_right_tail_prob_list <- function(mu_theta_hat_mat, max_needed, Omega) {
  # obtain the pmf list
  pmf_list <- apply(X = mu_theta_hat_mat, MARGIN = 1, FUN = function(curr_row) {
    # ensure the tails go out to at least 1e-16 but not farther than 1e-50.
    max_count <- max(max_needed, qnbinom(p = 1e-16, mu = curr_row[["mu"]], size = curr_row[["theta"]], lower.tail = FALSE))
    max_count <- min(max_count, qnbinom(p = 1e-50, mu = curr_row[["mu"]], size = curr_row[["theta"]], lower.tail = FALSE))
    dnbinom(x = seq(0L, max_count), mu = curr_row[["mu"]], size = curr_row[["theta"]])
  }, simplify = FALSE)

  # convert each omega to a binary string; compute its occupancy count; initialize the conv pmf list
  occupancy_counts <- rowSums(Omega)
  binary_str_labels <- apply(X = Omega, MARGIN = 1, FUN = function(r) paste0(r, collapse = ""))
  conv_pmf_list <- vector(mode = "list", length = nrow(Omega))
  names(conv_pmf_list) <- binary_str_labels

  # iterate over channels of different occupancy counts
  for (occupancy_count in sort(unique(occupancy_counts))) {
    omega_idxs <- which(occupancy_counts == occupancy_count)
    for (omega_idx in omega_idxs) {
      if (occupancy_count == 1L) {
        pmf_list_idx <- which(Omega[omega_idx,] == 1L)
        conv_pmf_list[[omega_idx]] <- pmf_list[[pmf_list_idx]]
      } else {
        binary_str_label <- binary_str_labels[omega_idx]
        binary_str_label_split <- strsplit(binary_str_label, split = "")[[1]]
        # find rightmost "1"
        rightmost_1_idx <- max(which(binary_str_label_split == "1"))
        # construct left and right convolution strings
        left_conv_string <- binary_str_label_split
        left_conv_string[rightmost_1_idx] <- "0"
        right_conv_string <- rep("0", length(left_conv_string))
        right_conv_string[rightmost_1_idx] <- "1"
        # index the left and right convolutions
        left_conv <- conv_pmf_list[[paste0(left_conv_string, collapse = "")]]
        right_conv <- conv_pmf_list[[paste0(right_conv_string, collapse = "")]]
        # convolve the left and right convolutions
        combined_conv <- convolve_pmfs_v2(a = left_conv, b = right_conv)
        conv_pmf_list[[binary_str_label]] <- combined_conv
      }
    }
  }

  # floor and normalize the convolved pmfs
  for (i in seq_along(conv_pmf_list)) {
    curr_pmf <- conv_pmf_list[[i]]
    curr_pmf <- pmax(curr_pmf, 1e-50)
    curr_pmf <- curr_pmf / sum(curr_pmf)
    conv_pmf_list[[i]] <- curr_pmf
  }

  # compute the right-tail probabilities
  right_tail_prob_list <- lapply(X = conv_pmf_list, FUN = function(curr_pmf) {
    rev(cumsum(rev(curr_pmf)))
  })
  return(right_tail_prob_list)
}
