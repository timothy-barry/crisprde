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
#' umi_tab_fp <- "/Users/timbarry/research_offsite/external/crispr-quant/guideseq/count_tables/293T-SpRY-Cas9-1620-GSPneg-S79_S87_L001_count_table.rds"
#' count_df <- readRDS(umi_tab_fp)
#' res <- find_guideseq_edit_sites_dm_sliding_window(count_df)
#'
#' umi_tab_fp <- "/Users/timbarry/research_offsite/external/crispr-quant/guideseq/count_tables/293T-SpRY-Cas9-1620-GSPneg-S79_S87_L001_count_table.rds"
#' count_df <- readRDS(umi_tab_fp)
#' res_trt <- find_guideseq_edit_sites(count_df)
find_guideseq_edit_sites_dm_sliding_window <- function(count_df, window_size = 25, multiplicity_adjustment = "BH",
                                                       multiplicity_alpha = 0.1, chrs_to_keep = seq(1L, 22L), fit = NULL,
                                                       plot_col = c("dodgerblue3", "firebrick")[1]) {
  count_df <- count_df |> dplyr::filter(chr %in% chrs_to_keep) |> dplyr::ungroup()
  n_loci <- nrow(count_df)
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
  ds <- lapply(X = unique(count_df$chr), FUN = function(curr_chr) {
    coords <- count_df |> dplyr::filter(chr == curr_chr) |> dplyr::pull(coord)
    if (length(coords) == 1) {
      1e7
    } else {
      dif_v <- diff(coords)
      first_d <- dif_v[1]
      last_d <- dif_v[length(dif_v)]
      middle_d <- pmin(dif_v[seq(1, length(dif_v) - 1)], dif_v[seq(2, length(dif_v))])
      c(first_d, middle_d, last_d)
    }
  }) |> unlist()
  count_df_w_dist$min_dist <- ds
  return(count_df_w_dist)
}
