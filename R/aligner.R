run_global_alignment <- function(candidate_protospacer, grna_spacer, match = 1, mismatch = 0, gopen = -1, gext = -1) {
  submat <- pwalign::nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly = TRUE)
  pairwiseAlignment(pattern = grna_spacer,
                    subject = candidate_protospacer,
                    substitutionMatrix = submat,
                    gapOpening = gopen,
                    gapExtension = gext,
                    type = "global",
                    scoreOnly = TRUE)
}

find_best_pam_on_strand <- function(query_seq, grna_spacer, pam_pattern) {
  spacer_length <- length(grna_spacer)
  pam_hits <- matchPattern(pam_pattern, query_seq, fixed = FALSE)
  if (length(pam_hits) == 0L) return(data.frame())
  pam_starts <- start(pam_hits)
  window_starts <- pam_starts - spacer_length
  ok_window <- window_starts >= 1L
  # subset hits to ok window
  pam_starts <- pam_starts[ok_window]; window_starts <- window_starts[ok_window]
  # obtain the subsetted query strings (i.e., candidate protospacers)
  views <- DNAStringSet(Views(query_seq, start = window_starts, width = spacer_length))
  scores <- vapply(views, run_global_alignment, numeric(1), grna_spacer = grna_spacer)
  data.frame(score = scores, pam_start = pam_starts)
}

align_spry_cas9 <- function(query_seq, grna_spacer, pam_pattern = "NGG") {
  plus_df <- find_best_pam_on_strand(query_seq = query_seq,
                                     grna_spacer = grna_spacer,
                                     pam_pattern = pam_pattern) |> dplyr::mutate(pam_strand = "+")
  plus_df$cut_base <- plus_df$pam_start - 3L
  minus_df <- find_best_pam_on_strand(query_seq = reverseComplement(query_seq),
                                      grna_spacer = grna_spacer,
                                      pam_pattern = pam_pattern) |> dplyr::mutate(pam_strand = "-")
  # translate minus-strand coordinates back to original system
  if (nrow(minus) >= 1) {
    minus_df$pam_start <- length(query_seq) - minus_df$pam_start - 1L
    minus_df$cut_base <- minus_df$pam_start + 5L
  }
  # construct combined df
  combined_df <- rbind(plus_df, minus_df)
}
