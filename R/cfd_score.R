calulate_cfd_score <- function(homology_dna, homology_gRNA, homology_has_hit) {
  mismatch_scores <- readRDS(system.file("extdata", "cfd_mismatch_scores.rds", package = "crisprde"))

  revcom <- function(s) {
    basecomp <- c(A = "T", C = "G", G = "C", T = "A", U = "A", "-" = "-")
    paste0(basecomp[rev(strsplit(s, "", fixed = TRUE)[[1]])], collapse = "")
  }

  score_one <- function(dna_seq, guide_seq, homology_has_hit) {
    if (!homology_has_hit) return(NA_real_)

    dna_seq <- toupper(dna_seq)
    guide_seq <- toupper(guide_seq)
    dna_protospacer <- substr(dna_seq, 1, nchar(dna_seq) - 3)
    dna_protospacer <- gsub("T", "U", dna_protospacer, fixed = TRUE)
    guide_seq <- gsub("T", "U", guide_seq, fixed = TRUE)
    score <- 1

    for (i in seq_len(nchar(dna_protospacer))) {
      guide_base <- substr(guide_seq, i, i)
      dna_base <- substr(dna_protospacer, i, i)

      if (guide_base != dna_base) {
        key <- paste0("r", guide_base, ":d", revcom(dna_base), ",", i)
        mismatch_score <- mismatch_scores[key]
        if (is.na(mismatch_score)) {
          mismatch_score <- 1
        }
        score <- score * mismatch_score
      }
    }
    return(unname(score))
  }
  unname(mapply(score_one, homology_dna, homology_gRNA, homology_has_hit, USE.NAMES = FALSE))
}
