run_global_alignment <- function(candidate_protospacer, grna_spacer, match = 1, mismatch = 0, gopen = -1, gext = -1) {
  submat <- pwalign::nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly = TRUE)
  pwalign::pairwiseAlignment(pattern = grna_spacer,
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


#' Align spacer sequence
#'
#' @param query_seq the target sequence (typically ~50-100 bp)
#' @param grna_spacer the 20-bp gRNA spacer sequence
#' @param pam_pattern the PAM (default: NGG)
#'
#' @returns a data frame summarizing the alignment results
#' @export
align_spacer_seq <- function(query_seq, grna_spacer, pam_pattern = "NGG") {
  grna_length <- length(grna_spacer)
  plus_df <- find_best_pam_on_strand(query_seq = query_seq,
                                     grna_spacer = grna_spacer,
                                     pam_pattern = pam_pattern) |> dplyr::mutate(pam_strand = "+")
  minus_df <- find_best_pam_on_strand(query_seq = reverseComplement(query_seq),
                                      grna_spacer = grna_spacer,
                                      pam_pattern = pam_pattern) |> dplyr::mutate(pam_strand = "-")
  if (nrow(plus_df) >= 1) {
    plus_df <- plus_df |> dplyr::mutate(pam_end = pam_start + 2L,
                                        cut_base_start = pam_start - 4L,
                                        cut_base_end = pam_start - 3L,
                                        protospacer_start = pam_start - grna_length,
                                        protospacer_end = pam_start - 1L)
  }
  # translate minus-strand coordinates back to original system
  if (nrow(minus_df) >= 1) {
    minus_df <- minus_df |> dplyr::mutate(pam_start = length(query_seq) - pam_start - 1L,
                                          pam_end = pam_start + 2L,
                                          cut_base_start = pam_start + 5L,
                                          cut_base_end = pam_start + 6L,
                                          protospacer_start = pam_start + 3L,
                                          protospacer_end = pam_start + 2L + grna_length)
  }
  # construct combined df
  combined_df <- rbind(plus_df, minus_df) |> dplyr::arrange(dplyr::desc(score))
}


#' Align a spacer sequence to the genome
#'
#' Aligns a short (20pb) spacer sequence to the genome.
#'
#' @param grna_spacer gRNA spacer sequence (written in 5' -> 3' direction)
#' @param pam (optional) the PAM sequence
#' @param ref_genome file path to the reference genome
#'
#' @returns a GAlignments object containing the alignmenet
#' @export
#'
#' @examples
#' alignment <- align_spacer_to_genome("ACTGATAGGGGTCGCGGTAG")
align_spacer_to_genome <- function(grna_spacer, ref_genome = .get_config_path("REF_GENOME"), pam = NULL) {
  if (is(grna_spacer, "DNAString")) {
    grna_spacer <- as.character(grna_spacer)
  }
  if (!is.null(pam)) {
    grna_spacer <- paste0(grna_spacer, pam)
  }

  temp_dir <- tempdir()
  temp_fasta <- paste0(temp_dir, "/spacer.fa")
  writeLines(c(">spacer", grna_spacer), temp_fasta)
  temp_sai <- paste0(temp_dir, "/spacer.sai")
  temp_sam <- paste0(temp_dir, "/spacer.sam")
  temp_bam <- paste0(temp_dir, "/spacer")
  bwa_cmd_1 <- paste0("bwa aln -n 0 -o 0 ", ref_genome, " ", temp_fasta," > ", temp_sai)
  bwa_cmd_2 <- paste0("bwa samse ", ref_genome, " ",  temp_sai, " ",  temp_fasta, " > ", temp_sam)
  bwa_cmd <- paste0(bwa_cmd_1, "; ", bwa_cmd_2)

  # run alignment
  system(bwa_cmd)

  # read SAM file
  Rsamtools::asBam(temp_sam, destination = temp_bam, overwrite = TRUE)
  alignment <- GenomicAlignments::readGAlignments(file = paste0(temp_bam, ".bam"))
  return(alignment)
}
