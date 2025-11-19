#' Find GUIDE-seq edit sites DM
#'
#' @param count_df
#' @param window_size
#' @param multiplicity_adjustment
#' @param chrs_to_keep
#' @param plot_col
#'
#' @returns
#' @export
#'
#' @examples
#' data_dir <- .get_config_path("LOCAL_CRISPR_DE_DATA_DIR")
#'
#' umi_tab_fp <- paste0(data_dir, "/guideseq/count_tables/293T-SpRY-Cas9-dsODN-only-GSPneg-S81_S89_L001_count_table.rds")
#' count_df <- readRDS(umi_tab_fp)
#' res_cntrl <- find_guideseq_edit_sites(count_df, bin_genome = FALSE, plot_col = "firebrick")
#'
#' umi_tab_fp <- paste0(data_dir, "guideseq/count_tables/293T-SpRY-Cas9-1620-GSPneg-S79_S87_L001_count_table.rds")
#' count_df <- readRDS(umi_tab_fp)
#' res_trt <- find_guideseq_edit_sites_dm_sliding_window(count_df)
find_guideseq_edit_sites_dm_sliding_window <- function(count_df, window_size = 25, multiplicity_adjustment = "BH",
                                                       multiplicity_alpha = 0.1, chrs_to_keep = seq(1L, 22L), fit = NULL,
                                                       plot_col = c("dodgerblue3", "firebrick")[1]) {
  count_df <- count_df |> dplyr::filter(chr %in% chrs_to_keep) |> dplyr::ungroup()
  n_loci <- nrow(count_df)

  if (TRUE) { # original
    # obtain X matrix
    X <- sapply(X = seq(1, n_loci), FUN = function(i) {
      curr_row <- count_df[i,]
      curr_coord <- curr_row$coord
      curr_chr <- curr_row$chr
      curr_count <- curr_row$count
      # subset count df to current chromosome
      offset <- (window_size - 3)/2
      count_df_chr_sub <- count_df |> dplyr::filter(chr == curr_chr)
      k_0 <- curr_count
      k_minus_1 <- count_df_chr_sub |> dplyr::filter(coord == (curr_coord - 1)) |> dplyr::pull(count) |> sum()
      k_plus_1 <- count_df_chr_sub |> dplyr::filter(coord == (curr_coord + 1)) |> dplyr::pull(count) |> sum()
      k_minus_2 <- count_df_chr_sub |> dplyr::filter(coord <= (curr_coord - 2),
                                                     coord >= (curr_coord - offset)) |> dplyr::pull(count) |> sum()
      k_plus_2 <- count_df_chr_sub |> dplyr::filter(coord >= (curr_coord + 2),
                                                    coord <= (curr_coord + offset)) |> dplyr::pull(count) |> sum()
      c(k_minus_2, k_minus_1, k_0, k_plus_1, k_plus_2)
    }) |> t()
  }

  # fit model
  if (is.null(fit)) {
    print("Fitting DM model.")
    fit <- fit_dm_model(X)
  }

  # compute p-values and update count_df
  p_values <- compute_lrt_p_vals(X = X, pi = fit$pi, theta = fit$theta)
  res_df <- count_df |> dplyr::mutate(p_value = p_values)

  # partition the genome into 100-bp windows; compute clustered result
  res_df <- compute_distances_between_occupied_loci(res_df, 20)
  clustered_res_df <- res_df |>
    dplyr::group_by(group_id) |>
    dplyr::summarize(p_combined = combine_p_values_simes(p_value))
  clustered_res_df$significant_hit <- p.adjust(p = clustered_res_df$p_combined, method = "BH") < 0.1

  # make plots
  # 1. manhattan plot
  manhattan_plot <- make_manhattan_plot(res_df)
  # 2. p-value histogram
  p_value_histogram <- make_p_value_histogram(clustered_res_df)
  # 3. global scatterplots
  global_scatterplot_linear <- make_scatterplot(count_df = count_df, x_range = NULL,
                                                facet_on_chr = TRUE, log_trans = FALSE, col = plot_col)
  global_scatterplot_log <- make_scatterplot(count_df = count_df, x_range = NULL,
                                             facet_on_chr = TRUE, log_trans = TRUE, col = plot_col)
  # 4. zoomed in plots of significant sites
  zoomed_plots <- make_discovery_site_scatterplots_dm(res_df, clustered_res_df, col = plot_col)

  out <- list(res_df = res_df, clustered_res_df = clustered_res_df,
              count_matrix = X, fit = fit, manhattan_plot = manhattan_plot,
              p_value_histogram = p_value_histogram,
              global_scatterplot_linear = global_scatterplot_linear,
              global_scatterplot_log = global_scatterplot_log,
              zoomed_plots = zoomed_plots)
  return(out)
}

#' Find GUIDE-seq edit sites DM clustered
#'
#' @param count_df
#' @param max_chain_link
#' @param multiplicity_adjustment
#' @param multiplicity_alpha
#' @param chrs_to_keep
#' @param fit
#' @param plot_col
#'
#' @returns
#' @export
#'
#' @examples
#' data_dir <- .get_config_path("LOCAL_CRISPR_DE_DATA_DIR")
#'
#' umi_tab_fp <- paste0(data_dir, "/guideseq/count_tables/293T-SpRY-Cas9-dsODN-only-GSPneg-S81_S89_L001_count_table.rds")
#' count_df <- readRDS(umi_tab_fp)
#' res_cntrl <- find_guideseq_edit_sites_dm_clustered(count_df, plot_col = "firebrick")
#'
#' umi_tab_fp <- paste0(data_dir, "guideseq/count_tables/293T-SpRY-Cas9-1620-GSPneg-S79_S87_L001_count_table.rds")
#' count_df <- readRDS(umi_tab_fp)
#' res_trt <- find_guideseq_edit_sites_dm_clustered(count_df)
find_guideseq_edit_sites_dm_clustered <- function(count_df, max_chain_link = 25, multiplicity_adjustment = "BH",
                                                  multiplicity_alpha = 0.1, chrs_to_keep = seq(1L, 22L), fit = NULL,
                                                  plot_col = c("dodgerblue3", "firebrick")[1]) {
  # 1. subset count_df, keeping only chromosomes in chrs_to_keep
  count_df <- count_df |> dplyr::filter(chr %in% chrs_to_keep)

  # 2. get the clustered count df
  distance_df <- compute_distances_between_occupied_loci(count_df = count_df, thresh = max_chain_link)

  # 3. generate X matrix
  X <- distance_df |>
    dplyr::group_by(group_id) |>
    dplyr::do(construct_x_row(.))
  X_mat <- as.matrix(X[,c("x1", "x2", "x3", "x4", "x5")])

  # 4. fit dm model
  if (is.null(fit)) {
    print("Fitting DM model.")
    fit <- fit_dm_model(X_mat)
    keep <- apply(X_mat, MARGIN = 1, FUN = function(x) sum(x >= 1)) >= 2
    fit_thin <- fit_dm_model(X_mat[keep,])
  }

  # 5. compute p-values and update df
  p_vals <- compute_lrt_p_vals(X_mat, pi = fit$pi, theta = fit$theta)
  X$p_value <- p_vals
  X$significant_hit <- p.adjust(p = p_vals, method = "BH") < multiplicity_alpha
  # for plotting
  X$coord <- X$lead_base
  X$p_combined <- X$p_value
  X <- X |> dplyr::arrange(p_value)

  # 6. make plots
  # a. manhattan plot
  manhattan_plot <- make_manhattan_plot(X)
  # b. p-value histogram
  p_value_histogram <- make_p_value_histogram(X)
  # c. global scatterplots
  global_scatterplot_linear <- make_scatterplot(count_df = count_df, x_range = NULL,
                                                facet_on_chr = TRUE, log_trans = FALSE, col = plot_col)
  global_scatterplot_log <- make_scatterplot(count_df = count_df, x_range = NULL,
                                             facet_on_chr = TRUE, log_trans = TRUE, col = plot_col)
  # d. zoomed in plots of significant sites
  zoomed_plots <- make_discovery_site_scatterplots_dm(distance_df, X, col = plot_col)
  X$coord <- X$p_combined <- NULL

  out <- list(distance_df = distance_df,
              X = X, fit = fit,
              manhattan_plot = manhattan_plot,
              p_value_histogram = p_value_histogram,
              global_scatterplot_linear = global_scatterplot_linear,
              global_scatterplot_log = global_scatterplot_log,
              zoomed_plots = zoomed_plots)
}

construct_x_row <- function(df, window_size = 25) {
  lead_base_pos <- which.max(df$count)
  lead_base <- df$coord[lead_base_pos]
  offset <- (window_size - 3)/2
  x1 <- df |> dplyr::filter(coord >= (lead_base - offset), coord <= (lead_base - 2)) |> dplyr::pull(count) |> sum()
  x2 <- df |> dplyr::filter(coord == (lead_base - 1)) |> dplyr::pull(count) |> sum()
  x3 <- df[[lead_base_pos, "count"]]
  x4 <- df |> dplyr::filter(coord == (lead_base + 1)) |> dplyr::pull(count) |> sum()
  x5 <- df |> dplyr::filter(coord <= (lead_base + offset), coord >= (lead_base + 2)) |> dplyr::pull(count) |> sum()
  with(df, data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5,
                      chr = chr[1], range_start = coord[1], range_end = coord[length(coord)],
                      lead_base = lead_base)) |>
    dplyr::mutate(region_str = paste0("chr", chr, ":", range_start, "-", range_end))
}

combine_p_values_simes <- function(ps, alpha = 0.1) {
  sorted_ps <- sort(ps)
  n <- length(sorted_ps)
  min(n * sorted_ps/seq(1, n))
}

# clustering loci; computing distances between occupied loci
compute_distances_between_occupied_loci <- function(count_df, thresh = 20) {
  # 1. simple distance
  count_df_w_dist <- count_df |>
    dplyr::group_by(chr) |>
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
  count_df_w_dist <- count_df_w_dist |>
    dplyr::ungroup() |> dplyr::mutate(group_id = group_id)

  # 2. min distance
  #ds <- lapply(X = unique(count_df$chr), FUN = function(curr_chr) {
  #  coords <- count_df |> dplyr::filter(chr == curr_chr) |> dplyr::pull(coord)
  #  if (length(coords) == 1) {
  #    1e7
  #  } else {
  #    dif_v <- diff(coords)
  #    first_d <- dif_v[1]
  #    last_d <- dif_v[length(dif_v)]
  #    middle_d <- pmin(dif_v[seq(1, length(dif_v) - 1)], dif_v[seq(2, length(dif_v))])
  #    c(first_d, middle_d, last_d)
  #  }
  #}) |> unlist()
  #count_df_w_dist$min_dist <- ds
  return(count_df_w_dist)
}
