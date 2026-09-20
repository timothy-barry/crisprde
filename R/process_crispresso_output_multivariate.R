#' Process CRISPRESSO allele table
#'
#' This function assumes that CRISPRESSO has been run with stringent settings, in particular with `--ignore_substitutions`, `quantification_window_size 1`, `--quantification_window_center -3`
#'
#' @param allele_table_fp the file path to the allele count table
#'
#' @returns a data frame with columns `allele_id` (giving the aligned sequence concatenated to the reference sequence) `mutation_type `(insertion, deletion, or complex), `mutation_length`, `read_count`, and `mutated`. Alleles classified by CRISPResso as "unmodified" are collapsed into a single allele labeled "unmutated."
#' @export
#'
#' @examples
#' allele_table_fp <- "/Users/timbarry/research_offsite/external/bauer-lab/rhampseq_bcl11a/hiseq1_crispresso_output/GE0423-edited/CRISPResso_on_1450_OT_0000/Alleles_frequency_table.txt"
#' allele_feature_table <- process_crispresso_allele_table_cas9(allele_table_fp)
#'
#' # TO-DO: try to add position of the mutation relative to the cut site (whatever that might mean).
#' # Note: this field would require us to know the cutsite.
process_crispresso_allele_table_cas9 <- function(allele_table_fp) {
  allele_tab <- readr::read_delim(allele_table_fp) |>
    dplyr::select(aligned_seq = Aligned_Sequence, ref_seq = Reference_Sequence,
                  status = Read_Status, n_deleted, n_inserted, n_reads = "#Reads")
  allele_tab$allele_id <- paste0(allele_tab$aligned_seq, ":", allele_tab$ref_seq)

  # 1. separate mutated from unmutated alleles (as determined by CRISPResso)
  unmutated_tab <- allele_tab |> dplyr::filter(status == "UNMODIFIED")
  mutated_tab <- allele_tab |> dplyr::filter(status == "MODIFIED")

  if (nrow(mutated_tab) >= 1L) {
    # 3. determine mutation type
    mutation_types_new <- sapply(X = seq_len(nrow(mutated_tab)), FUN = function(i) {
      deletion_present <- mutated_tab$n_deleted[i] >= 1
      insertion_present <- mutated_tab$n_inserted[i] >= 1
      if (deletion_present & !insertion_present) {
        mutation_type <- "deletion"
      } else if (!deletion_present & insertion_present) {
        mutation_type <- "insertion"
      } else {
        mutation_type <- "complex"
      }
    })

    # 4. determine mutation length
    mutation_lengths_new <- mutated_tab$n_deleted + mutated_tab$n_inserted
    mutated_tab_out <- data.frame(allele_id = mutated_tab$allele_id,
                                  mutation_type = mutation_types_new, mutation_length = mutation_lengths_new,
                                  read_count = mutated_tab$n_reads, mutated = TRUE)
  } else {
    mutated_tab_out <- data.frame()
  }

  # 6. prepare the output data frame
  out_df <- rbind(
    data.frame(allele_id = "unmutated", mutation_type = NA_character_,
               mutation_length = NA_integer_,
               read_count = sum(unmutated_tab$n_reads),
               mutated = FALSE),
    mutated_tab_out
    )
  return(out_df)
}
