
boost_p_values_ihw <- function(augmented_result_df, multiplicity_alpha = 0.2) {
  # compute the align score as homology_n_mismatches + 2 * homology_n_bulges; for na windows, the max of this quantity + 1
  augmented_result_df <- augmented_result_df |>
    dplyr::mutate(align_score = homology_n_mismatches + 2 * homology_n_bulges)
  max_align_score <- augmented_result_df$align_score |> max(na.rm = TRUE)
  augmented_result_df$align_score[is.na(augmented_result_df$align_score)] <- max_align_score + 1
  ihw_fit <- IHW::ihw(pvalues = augmented_result_df$p_value, covariates = -augmented_result_df$align_score,
                      alpha = multiplicity_alpha, nbins = 4)
  w <- IHW::weights(ihw_fit)
  augmented_result_df$p_weighted <- augmented_result_df$p_value/w
  return(augmented_result_df)
}


#' Boost p-values (Genovese)
#'
#' @param augmented_result_df output of `overlap_homology_and_result_dfs`
#' @param multiplicity_alpha nominal FDR
#' @param gamma gamma
#'
#' @examples
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' count_df <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#'  dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & treated & replicate_id %in% 1:2) |>
#'  dplyr::filter(chr != "chrM") |>
#'  dplyr::select(chr, coord, strand, umi_count, primer_type, replicate_id)
#' Y_mat <- construct_replicate_count_table(count_df)
#' result_df <- run_multireplicate_guideseq_method(Y_mat = Y_mat, lambda = 10, c_tukey_sigma = 50, multiplicity_alpha = 0.2, robust_fit = TRUE, incorporate_occupancy_info = TRUE)$res_df
#' homology_df <- load_crispritz_output("/Users/timbarry/research_offsite/external/bauer-lab/guideseq_elane/crispritz/crispritz_CCCCGGCAGAAACGTCCGCG.hg38.targets.txt")
#' augmented_result_df <- overlap_homology_and_result_dfs(result_df, homology_df)
#' weighted_result_df <- boost_p_values_genovese(augmented_result_df)
#' qq_plot <- weighted_result_df |> dplyr::mutate(p_value = p_value_weighted) |> make_guideseq_qq_plot()
boost_p_values_genovese <- function(augmented_result_df, multiplicity_alpha = 0.2, gamma = 0.5) {
  # compute the alignment score
  align_score <- augmented_result_df$homology_n_mismatches + 2 * augmented_result_df$homology_n_bulges
  align_score[is.na(align_score)] <- max(align_score, na.rm = TRUE) + 1
  # compute weights
  w <- exp(-gamma * align_score)
  w_tilde <- w/mean(w)
  p_value_weighted <- augmented_result_df$p_value/w_tilde
  q_value_weighted <- p.adjust(p = p_value_weighted, method = "BH")
  nominated_window_weighted <- q_value_weighted < multiplicity_alpha
  out <- augmented_result_df |> dplyr::mutate(p_value_weighted = pmin(1, p_value_weighted),
                                       nominated_window_weighted = nominated_window_weighted) |>
    dplyr::arrange(p_value_weighted)
  return(out)
}
