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
#' umi_tab_fp <- "/Users/timbarry/research_offsite/external/crispr-quant/guideseq/count_tables/293T-SpRY-Cas9-dsODN-only-GSPneg-S81_S89_L001_count_table.rds"
#' count_df <- readRDS(umi_tab_fp)
find_guideseq_edit_sites_dm_sliding_window <- function(count_df, window_size = 25, multiplicity_adjustment = "BH",
                                                       chrs_to_keep = seq(1L, 22L), plot_col = c("dodgerblue3", "firebrick")[1]) {
  count_df <- count_df |> dplyr::filter(chr %in% chrs_to_keep)
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
  fit <- fit_dm_model(X)
  # compute p-values
  p_values <- compute_lrt_p_vals(X = X, pi = fit$pi, theta = fit$theta)
  count_df$p_value <- p_values
  out <- list(count_df = count_df, count_matrix = X, fit = fit)
  return(out)
}
