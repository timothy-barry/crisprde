#' Tune hyperparameters
#'
#' @param Y_mat_trt treated count matrix
#' @param Y_mat_cntrl control count matrix
#' @param c_grid robust hyperparm grid
#' @param lambda_grid lambda grid
#' @param incorporate_occupancy_info a boolean (T/F) indicating whether to incorporate occupancy information into the p-value calculation
#' @param count_model Either "nb" or "nb_mixture" for all channels, or a named
#'   character vector assigning one of these models to each channel (matrix row name).
#'   The same mapping applies to both conditions. Both models use robust estimation
#'   and share the tuning grid for c.
#' @param tail_start Upper bound on expected null p-values used for the control tail KS statistic.
#' @param ks_tolerance Maximum acceptable tail KS as a multiple of the minimum across the grid.
#' @param weight_p_values Whether to add `p_value_weighted` and `p_value_weight` to the results; `p_value` always contains raw p-values.
#' @param annotated_clustered_count_df_trt annotated clustered count data frame for the treated condition
#' @param annotated_clustered_count_df_cntrl annotated clustered count data frame for the control condition
#' @param tau baseline normalized weight for zero homology scores
#' @param verbose whether to print detailed progress messages within steps
#'
#' @returns A list containing `result_list` (skinny window-level tables indexed by
#'   condition, c, and lambda), `nb_model_fits` (parameter matrices indexed by
#'   condition and c, with unused parameter columns set to NA for mixed models),
#'   `occupancy_model_fits` (occupancy fits indexed by `trt`
#'   and `cntrl`, including `pi_hat` and `tbp_pattern_df`),
#'   `selected_params`, `selected_trt_run`, and `selected_cntrl_run`. Selected
#'   runs contain full annotations sorted by raw `p_value` and their fitted count
#'   parameters. When weighting is enabled, results also contain `p_value_weighted`
#'   and `p_value_weight`. Selection uses raw p-values for control calibration
#'   and treated on-target sensitivity; no sites are nominated.
#'   Occupancy fits are returned independently of selection. The returned
#'   `count_model` retains the supplied model label or channel-to-model mapping.
#' @export
#'
#' @examples
#' # with p-value boosting and filtering on homology
#' #######
#' # ELANE
#' #######
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' count_df_all <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#' dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & replicate_id %in% 1:2, chr != "chrM")
#' homology_df <- load_crispritz_output("/Users/timbarry/research_offsite/external/bauer-lab/guideseq_elane/crispritz_CCCCGGCAGAAACGTCCGCG.hg38.targets.txt")
#' on_target_posit <- readRDS(paste0(elane_dir, "spacer_info.rds"))$elane_posit
#' n_run_df <- load_n_run_bed("/Users/timbarry/research_offsite/ref_genome_dir/hg38_N_runs_min10.bed")
#' annotated_clustered_count_df_trt <- count_df_all |> dplyr::filter(treated) |> cluster_loci() |>
#'   annotate_clustered_count_df(homology_df = homology_df, n_run_df = n_run_df, on_target_posit = on_target_posit)
#' annotated_clustered_count_df_cntrl <- count_df_all |> dplyr::filter(!treated) |> cluster_loci() |>
#'   annotate_clustered_count_df(homology_df = homology_df, n_run_df = n_run_df)
#' Y_mat_trt <- construct_replicate_count_table(annotated_clustered_count_df_trt)
#' Y_mat_cntrl <- construct_replicate_count_table(annotated_clustered_count_df_cntrl)
#' hyperparam_res <- tune_hyperparameters(Y_mat_trt = Y_mat_trt, Y_mat_cntrl = Y_mat_cntrl,
#'   annotated_clustered_count_df_trt = annotated_clustered_count_df_trt,
#'   annotated_clustered_count_df_cntrl = annotated_clustered_count_df_cntrl)
#'
#' ################
#' # SBDSP1 DONOR 1
#' ################
#' options(parallelly.fork.enable = TRUE)
#' future::plan(future::multicore, workers = 5)
#' sbds_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_sbds_ns_15/")
#' count_df_all <- readRDS(paste0(sbds_dir, "count_tables_without_multimap/combined_count_df.rds")) |>
#'   dplyr::filter(replicate_id == 1L, chr != "chrM")
#' homology_df <- load_crispritz_output(paste0(sbds_dir, "crispritz_ACTGAAATCTGTAAGCAGGC_hg38.targets.txt"))
#' on_target_posit <- readRDS(paste0(sbds_dir, "spacer_info.rds"))$sbdsp1_posit
#' n_run_df <- load_n_run_bed("/Users/timbarry/research_offsite/ref_genome_dir/hg38_N_runs_min10.bed")
#' blacklist_df <- load_encode_blacklist_bed("/Users/timbarry/research_offsite/ref_genome_dir/hg38-blacklist.v2.bed")
#' annotated_clustered_count_df_trt <- count_df_all |> dplyr::filter(treated) |> cluster_loci() |>
#'   annotate_clustered_count_df(homology_df = homology_df, n_run_df = n_run_df,
#'                              encode_blacklist_df = blacklist_df,
#'                              on_target_posit = on_target_posit) |>
#'   dplyr::filter(!overlaps_n_run, !overlaps_encode_blacklist)
#' annotated_clustered_count_df_cntrl <- count_df_all |> dplyr::filter(!treated) |> cluster_loci() |>
#'   annotate_clustered_count_df(homology_df = homology_df, n_run_df = n_run_df,
#'                              encode_blacklist_df = blacklist_df) |>
#'   dplyr::filter(!overlaps_n_run, !overlaps_encode_blacklist)
#' Y_mat_trt <- construct_replicate_count_table(annotated_clustered_count_df_trt, channel_axes = "primer_type")
#' Y_mat_cntrl <- construct_replicate_count_table(annotated_clustered_count_df_cntrl, channel_axes = "primer_type")
#' hyperparam_res <- tune_hyperparameters(Y_mat_trt = Y_mat_trt, Y_mat_cntrl = Y_mat_cntrl,
#'   annotated_clustered_count_df_trt = annotated_clustered_count_df_trt,
#'   annotated_clustered_count_df_cntrl = annotated_clustered_count_df_cntrl,
#'   count_model = "nb_mixture")
#' fitted_count_dist_plot <- make_fitted_count_distribution_plots(
#'   Y_mat_trt = Y_mat_trt, Y_mat_cntrl = Y_mat_cntrl,
#'   hyperparam_res = hyperparam_res, c = 50)
tune_hyperparameters <- function(Y_mat_trt, Y_mat_cntrl,
                                 annotated_clustered_count_df_trt,
                                 annotated_clustered_count_df_cntrl,
                                 c_grid = c(5, 10, 25, 50, 100, 500, 1000),
                                 lambda_grid = c(0, 10, 25, 50, 100),
                                 count_model = c("nb", "nb_mixture")[1],
                                 incorporate_occupancy_info = TRUE,
                                 tail_start = 0.05, ks_tolerance = 1.5,
                                 weight_p_values = TRUE,
                                 lambda_default = 20, tau = 0.1,
                                 verbose = FALSE) {
  ############
  # CHECK ARGS
  ############
  if (!all(count_model %in% c("nb", "nb_mixture"))) {
    stop("`count_model` must be `nb` or `nb_mixture`.")
  }

  ###########################################################
  # PART 1: FIT OCCUPANCY AND COUNT MODELS TO BOTH CONDITIONS
  ###########################################################
  message("Step 1/5: Fitting occupancy and count models")
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
    if (verbose) message("Fitting count models for ", curr_condition)
    Y_mat <- Y_mat_list[[curr_condition]]
    channel_models <- get_channel_count_models(count_model, rownames(Y_mat))
    fits_by_model <- lapply(unique(channel_models), function(model) {
      fit_function <- switch(model,
        nb = fit_multirep_guideseq_count_null,
        nb_mixture = fit_multirep_guideseq_mixture_count_null)
      fit_function(Y_mat = Y_mat[channel_models == model, , drop = FALSE], c_grid = c_grid)
    })
    lapply(as.character(c_grid), function(c) {
      parameter_mats <- lapply(fits_by_model, function(fits) fits[[c]])
      parameters <- dplyr::bind_rows(lapply(parameter_mats, as.data.frame)) |> as.matrix()
      rownames(parameters) <- unlist(lapply(parameter_mats, rownames), use.names = FALSE)
      parameters[rownames(Y_mat), , drop = FALSE]
    }) |> setNames(c_grid)
  }) |> setNames(condition_grid)

  ##############################################################
  # PART 2: GROUP OBSERVATIONS AND COMPUTE TEST STATISTICS
  ##############################################################
  message("Step 2/5: Grouping observations and computing test statistics")
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
  message("Step 3/5: Computing p-values")
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
                                                       Omega = occupancy_fit$Omega,
                                                       count_model = count_model)
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

  ################################
  # PART 4: PREPARE WINDOW RESULTS
  ################################
  message("Step 4/5: Preparing window results")
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
      dplyr::select(window, dplyr::starts_with(c("homology", "window", "overlaps")),
                    if (condition == "trt") "on_target" else NULL) |>
      dplyr::distinct()
    result_df <- dplyr::left_join(x = result_df, y = annotation_df, by = "window")
  }) |> setNames(condition_grid)

  # expand p-values to windows and apply homology weighting
  result_list <- lapply(X = condition_grid, FUN = function(condition) {
    result_df_skinny <- result_dfs[[condition]] |> dplyr::select(window, group_id, homology_alignment_score,
                                                                 if (condition == "trt") "on_target" else NULL)
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
                                                          tau = tau)
        }
        result_df_skinny
      }) |> setNames(lambda_grid)
    }) |> setNames(c_grid)
  }) |> setNames(condition_grid)

  ##############################
  # PART 5: SUMMARIZE AND SELECT
  ##############################
  message("Step 5/5: Summarizing results and selecting hyperparameters")
  selected_params <- select_c_lambda_combination(result_list = result_list, c_grid = c_grid, lambda_grid = lambda_grid,
                                                 tail_start = tail_start, ks_tolerance = ks_tolerance)

  # join metadata only for the selected runs
  selected_c <- as.character(selected_params$c)
  selected_lambda <- as.character(selected_params$lambda)
  selected_runs <- lapply(X = condition_grid, FUN = function(condition) {
    res_df <- result_list[[condition]][[selected_c]][[selected_lambda]]
    ests_list <- list(mu_theta_hat_mat = nb_model_fits[[condition]][[selected_c]])
    res_df$lambda <- if (use_occupancy[[condition]]) selected_params$lambda else NA
    metadata_df <- result_dfs[[condition]] |>
      dplyr::select(-group_id, -homology_alignment_score, -dplyr::any_of("on_target"))
    res_df <- res_df |>
      dplyr::select(-group_id) |>
      dplyr::left_join(metadata_df, by = "window") |>
      dplyr::relocate(window, p_value, umi_count, lambda, occupancy_pattern) |>
      dplyr::arrange(p_value)
    list(res_df = res_df, ests_list = ests_list)
  }) |> setNames(condition_grid)
  ret <- list(selected_params = selected_params,
              selected_trt_run = selected_runs$trt,
              selected_cntrl_run = selected_runs$cntrl,
              result_list = result_list, nb_model_fits = nb_model_fits,
              occupancy_model_fits = occupancy_fit_list,
              count_model = count_model)

  return(ret)
}


select_c_lambda_combination <- function(result_list, c_grid, lambda_grid, tail_start = 0.05, ks_tolerance = 1.5) {
  # compute ks stat over cntrl dfs and on-target p-value for trt dfs
  result_list_trt <- result_list$trt
  result_list_cntrl <- result_list$cntrl

  # set up expected p-values
  n <- nrow(result_list_cntrl[[1]][[1]])
  expected_p <- ppoints(n = n)
  idxs_to_keep <- (expected_p < tail_start)
  expected_p_tail <- expected_p[idxs_to_keep]

  # compute the control-condition ks stat over the grid
  cntrl_tail_ks_stat <- lapply(X = c_grid, FUN = function(c) {
    lapply(X = lambda_grid, FUN = function(lambda) {
      res_df <- result_list_cntrl[[as.character(c)]][[as.character(lambda)]]
      observed_p <- sort(res_df$p_value)
      observed_p_tail <- observed_p[idxs_to_keep]
      ks_tail <- max(abs(observed_p_tail - expected_p_tail))
      c(c = c, lambda = lambda, ks_tail = ks_tail)
    }) |> dplyr::bind_rows()
  }) |> dplyr::bind_rows()

  # compute the treated-condition on-target p-value over the grid
  trt_on_target_p_val <- lapply(X = c_grid, FUN = function(c) {
    lapply(X = lambda_grid, FUN = function(lambda) {
      res_df <- result_list_trt[[as.character(c)]][[as.character(lambda)]]
      on_target_p <- res_df |> dplyr::filter(on_target) |> dplyr::pull(p_value)
      c(c = c, lambda = lambda, on_target_p = on_target_p)
    }) |> dplyr::bind_rows()
  }) |> dplyr::bind_rows()
  metric_df <- dplyr::left_join(cntrl_tail_ks_stat, y = trt_on_target_p_val, by = c("c", "lambda"))
  min_ks_tail <- min(metric_df$ks_tail)
  metric_df$acceptable_config <- (metric_df$ks_tail <= ks_tolerance * min_ks_tail)
  chosen_params <- metric_df |> dplyr::filter(acceptable_config) |> dplyr::arrange(on_target_p) |> dplyr::slice(1)
  ret <- data.frame(c = chosen_params$c, lambda = chosen_params$lambda)
  return(ret)
}


get_channel_count_models <- function(count_model, channels) {
  if (length(count_model) == 1L) {
    setNames(rep(count_model, length(channels)), channels)
  } else {
    count_model[channels]
  }
}


get_right_tail_prob_list <- function(mu_theta_hat_mat, max_needed, Omega, count_model) {
  # obtain the pmf list in channel order
  channel_models <- get_channel_count_models(count_model, rownames(mu_theta_hat_mat))
  pmf_list <- lapply(seq_len(nrow(mu_theta_hat_mat)), function(i) {
    curr_row <- mu_theta_hat_mat[i, ]
    if (channel_models[[i]] == "nb") {
      # clamp max_needed to fall between the 1e-16 and 1e-50 quantiles
      quantile_1eminus16 <- qnbinom(p = 1e-16, mu = curr_row[["mu"]], size = curr_row[["theta"]], lower.tail = FALSE)
      quantile_1eminus50 <- qnbinom(p = 1e-50, mu = curr_row[["mu"]], size = curr_row[["theta"]], lower.tail = FALSE)
      max_count <- min(max(max_needed, quantile_1eminus16), quantile_1eminus50)
      dnbinom(x = seq(0L, max_count), mu = curr_row[["mu"]], size = curr_row[["theta"]])
    } else {
      quantiles_1eminus16_1eminus50 <- sapply(X = c(1e-16, 1e-50), FUN = function(curr_quantile) {
        component_quantiles <- qnbinom(p = curr_quantile,
                                       mu = curr_row[c("mu_lower", "mu_upper")],
                                       size = curr_row[c("theta_lower", "theta_upper")],
                                       lower.tail = FALSE) |> max()
      })
      max_count <- min(max(max_needed, quantiles_1eminus16_1eminus50[1]), quantiles_1eminus16_1eminus50[2])
      lower_pmf <- dnbinom(x = seq(0L, max_count),
                           mu = curr_row[["mu_lower"]], size = curr_row[["theta_lower"]])
      upper_pmf <- dnbinom(x = seq(0L, max_count),
                           mu = curr_row[["mu_upper"]], size = curr_row[["theta_upper"]])
      pi <- curr_row[["pi"]]
      mixture_pmf <- (1 - pi) * lower_pmf + pi * upper_pmf
      return(mixture_pmf)
    }
  })

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
