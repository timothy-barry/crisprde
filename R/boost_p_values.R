#' Boost p-values using continuous homology weights
#'
#' @param augmented_result_df result data frame with homology annotations
#' @param multiplicity_alpha nominal FDR
#' @param tau baseline normalized weight for zero homology scores, between zero and one
#' @param gamma exponential distance-decay coefficient; use the same value for annotation
#'
#' @details Alignment scores are CFD * exp(-gamma * distance). Normalized weights
#'   are tau + A * score, where A = n * (1 - tau) / sum(score) and n is the
#'   number of windows. If all scores are zero, weights are one. Missing CFD
#'   or distance values contribute zero alignment score.
boost_p_values_genovese_cfd <- function(augmented_result_df, multiplicity_alpha = 0.5,
                                       tau = 0.1, gamma = log(20)/7) {
  w <- compute_alignment_scores(cfds = augmented_result_df$homology_cfd,
                                distances = augmented_result_df$homology_modal_base_cut_distance,
                                gamma = gamma)
  n <- length(w)
  if (sum(w) > 0) {
    A <- n * (1 - tau)/sum(w)
    w_tilde <- tau + A * w
  } else {
    w_tilde <- rep(1, n)
  }

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

compute_alignment_scores <- function(cfds, distances, gamma = log(20)/7) {
  w <- cfds * exp(-gamma * distances)
  w[is.na(w)] <- 0
  return(w)
}
