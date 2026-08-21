#' Boost p-values (Genovese)
#'
#' @param augmented_result_df result data frame with homology annotations
#' @param multiplicity_alpha nominal FDR
#' @param gamma_align tuning parameter controlling the homology-score penalty
#' @param gamma_distance tuning parameter controlling the modal-base-distance penalty
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
boost_p_values_genovese <- function(augmented_result_df, multiplicity_alpha = 0.5, gamma_align = NULL, gamma_distance = NULL) {
  # get the alignment score
  MAX_ALIGN_SCORE <- 9L
  align_score <- augmented_result_df$homology_n_mismatches + 2L * augmented_result_df$homology_n_bulges
  align_score <- pmin(align_score, MAX_ALIGN_SCORE)
  align_score[is.na(align_score)] <- MAX_ALIGN_SCORE + 1L
  # compute the distance of modal base to cut site
  MAX_MODAL_DISTANCE <- 19L
  modal_distance <- pmin(MAX_MODAL_DISTANCE, augmented_result_df$homology_modal_base_cut_distance)
  modal_distance[is.na(modal_distance)] <- MAX_MODAL_DISTANCE + 1L

  # compute weights
  if (is.null(gamma_align)) gamma_align <- log(50)/(2 * max(align_score))
  if (is.null(gamma_distance)) gamma_distance <- log(50)/(2 * max(modal_distance))
  w <- exp(-gamma_align * align_score - gamma_distance * modal_distance)
  w_tilde <- w/mean(w)

  # compute weighted p-values and discovery set
  p_value_weighted <- augmented_result_df$p_value/w_tilde
  q_value_weighted <- p.adjust(p = p_value_weighted, method = "BH")
  nominated_window_weighted <- q_value_weighted < multiplicity_alpha

  out <- augmented_result_df |>
    dplyr::mutate(p_value_unweighted = p_value,
                  nominated_window_unweighted = nominated_window,
                  p_value = pmin(1, p_value_weighted),
                  p_value_weight = w_tilde,
                  nominated_window = nominated_window_weighted) |>
    dplyr::arrange(p_value_weighted)
  return(out)
}


boost_p_values_genovese_cfd <- function(augmented_result_df, multiplicity_alpha = 0.5, prior_strength = "aggressive") {
  w <- compute_alignment_scores(cfds = augmented_result_df$homology_cfd,
                           distances = augmented_result_df$homology_modal_base_cut_distance,
                           prior_strength = prior_strength)
  prior_weights <- get_prior_weights(prior_strength = prior_strength)
  w[is.na(w)] <- min(prior_weights$cfd_weights) * min(prior_weights$distance_weights)
  w_tilde <- w/mean(w)

  # compute weighted p-values and discovery set
  p_value_weighted <- augmented_result_df$p_value/w_tilde
  q_value_weighted <- p.adjust(p = p_value_weighted, method = "BH")
  nominated_window_weighted <- q_value_weighted < multiplicity_alpha

  out <- augmented_result_df |>
    dplyr::mutate(p_value_unweighted = p_value,
                  nominated_window_unweighted = nominated_window,
                  p_value = pmin(1, p_value_weighted),
                  p_value_weight = w_tilde,
                  nominated_window = nominated_window_weighted) |>
    dplyr::arrange(p_value)
}

compute_alignment_scores <- function(cfds, distances, prior_strength = "aggressive") {
  prior_weights <- get_prior_weights(prior_strength = prior_strength)
  cfd_weights <- prior_weights$cfd_weights
  distance_weights <- prior_weights$distance_weights

  # cfd weighting
  cfd_breaks <- c(1, 0.8, 0.4, 0.2, 0.05, 0.01, -1)
  cfd_cuts <- cut(x = cfds, breaks = cfd_breaks)
  my_cfd_weights <- cfd_weights[as.integer(cfd_cuts)]

  # distance weighting
  distance_breaks <- c(-1, 0, 1, 2, 4, 6, Inf)
  distance_cuts <- cut(distances, breaks = distance_breaks)
  my_distance_weights <- distance_weights[as.integer(distance_cuts)]

  # final (normalized) weights
  w <- my_cfd_weights * my_distance_weights
  return(w)
}

get_prior_weights <- function(prior_strength = "aggressive") {
  if (prior_strength == "aggressive") {
    cfd_weights <- c(0.01, 0.05, 0.1, 0.25, 0.5, 1)
    distance_weights <- c(1, 0.5, 0.3, 0.1, 0.05, 0.01)
  } else if (prior_strength == "moderate") {
    cfd_weights <- c(0.05, 0.1, 0.2, 0.5, 0.8, 1)
    distance_weights <- c(1, 0.75, 0.5, 0.3, 0.1, 0.05)
  } else if (prior_strength == "mild") {
    cfd_weights <- c(0.25, 0.4, 0.6, 0.8, 0.9, 1)
    distance_weights <- c(1, 0.9, 0.8, 0.6, 0.4, 0.25)
  } else {
    stop("`prior_strength` not recognized. Choose one of `mild`, `moderate`, `aggressive`.")
  }
  out <- list(cfd_weights = cfd_weights, distance_weights = distance_weights)
}
