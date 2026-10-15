#' Boost p-values using continuous homology weights
#'
#' @param augmented_result_df result data frame containing `p_value` and `homology_alignment_score`
#' @param tau baseline normalized weight for zero homology scores, between zero and one
#'
#' @returns The input data frame with `p_value` unchanged, normalized weights
#'   in `p_value_weight`, and weighted p-values in `p_value_weighted`
#'   (the raw p-value divided by its weight, capped at one).
#'
#' @details Uses the stored homology_alignment_score. Normalized weights
#'   are tau + A * score, where A = n * (1 - tau) / sum(score) and n is the
#'   number of windows. If all scores are zero, weights are one. Missing scores
#'   contribute zero alignment score. No multiplicity correction or site nomination
#'   is performed.
boost_p_values_genovese_cfd <- function(augmented_result_df, tau = 0.1) {
  w <- augmented_result_df$homology_alignment_score
  w[is.na(w)] <- 0
  n <- length(w)
  if (sum(w) > 0) {
    A <- n * (1 - tau)/sum(w)
    w_tilde <- tau + A * w
  } else {
    w_tilde <- rep(1, n)
  }

  # compute weighted p-values
  p_value_weighted <- augmented_result_df$p_value/w_tilde

  out <- augmented_result_df |>
    dplyr::mutate(p_value_weighted = pmin(1, .env$p_value_weighted),
                  p_value_weight = w_tilde)
}

compute_alignment_scores <- function(cfds, distances, gamma = log(20)/7) {
  w <- cfds * exp(-gamma * distances)
  w[is.na(w)] <- 0
  return(w)
}
