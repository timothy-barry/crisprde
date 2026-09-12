#' Boost p-values (penalizing low-homology sites)
#'
#' @param augmented_result_df
#'
#' @param multiplicity_alpha
#' @param prior_strength
#'
#' @examples
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' clustered_count_df <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#'  dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & treated & replicate_id %in% 1:2 & chr != "chrM") |>
#'  dplyr::select(chr, coord, strand, umi_count, replicate_id) |>
#'  cluster_loci()
#' homology_df <- load_crispritz_output("/Users/timbarry/research_offsite/external/bauer-lab/guideseq_elane/crispritz_CCCCGGCAGAAACGTCCGCG.hg38.targets.txt")
#' annotated_clustered_count_df <- annotate_clustered_count_df(clustered_count_df, homology_df) # |> dplyr::filter(homology_has_hit)
#' Y_mat <- construct_replicate_count_table(annotated_clustered_count_df)
#' augmented_result_df <- run_multireplicate_guideseq_method(Y_mat = Y_mat, lambda = 10, c_tukey_sigma = 50, multiplicity_alpha = 0.2, robust_fit = TRUE, incorporate_occupancy_info = TRUE, annotated_clustered_count_df = annotated_clustered_count_df)$res_df
#'
#' # cfd weighting
#' weighted_result_df <- boost_p_values_genovese_cfd(augmented_result_df)
#' qq_plot <- weighted_result_df |> make_guideseq_qq_plot()
boost_p_values_genovese_cfd <- function(augmented_result_df, tau = 0.02, multiplicity_alpha = 0.5) {
  cfd_thresh <- 0.001
  distance_thresh <- 8L
  w_tilde <- numeric(length = nrow(augmented_result_df))

  # first, assign the low-homology penalty to the low-homology sites
  augmented_result_df <- augmented_result_df |>
    dplyr::mutate(low_homology = !homology_has_hit | homology_cfd < cfd_thresh | homology_modal_base_cut_distance >= distance_thresh)

  # derive the multiplicative constant for the high-homology sites
  n_low_homology_sites <- sum(augmented_result_df$low_homology)
  mult_constant <- nrow(augmented_result_df) - n_low_homology_sites * tau

  # compute the alignment weights among the high-homology windows
  w <- augmented_result_df |> dplyr::filter(!low_homology) |> dplyr::pull(homology_alignment_score)
  w_tilde_high_homology <- w * mult_constant / sum(w)

  # attach weights to augmented_result_df
  w_tilde <- numeric(length = nrow(augmented_result_df))
  w_tilde[augmented_result_df$low_homology] <- tau
  w_tilde[!augmented_result_df$low_homology] <- w_tilde_high_homology

  # apply BH
  p_value_weighted <- augmented_result_df$p_value/w_tilde
  q_value_weighted <- p.adjust(p = p_value_weighted, method = "BH")
  nominated_window_weighted <- q_value_weighted < multiplicity_alpha

  # prepare output
  out <- augmented_result_df |>
    dplyr::mutate(p_value_unweighted = p_value,
                  nominated_window_unweighted = nominated_window,
                  p_value = pmin(1, p_value_weighted),
                  p_value_weight = w_tilde,
                  nominated_window = nominated_window_weighted) |>
    dplyr::arrange(p_value)
}


compute_alignment_scores <- function(cfds, distances) {
  distance_thresh <- 8L
  max_distance_rel_weight <- 0.01
  coef <- log(max_distance_rel_weight)/(distance_thresh - 1)
  w <- exp(coef * distances) * cfds
  return(w)
}
