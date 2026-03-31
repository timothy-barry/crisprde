#' Process CRISPRESSO allele table
#'
#' This function assumes that CRISPRESSO has been run with stringent settings, in particular with `--ignore_substitutions`, `quantification_window_size 1`, `--quantification_window_center -3`
#'
#' @param allele_table_fp the file path to the allele count table
#'
#' @returns a data frame with columns `allele_id` (giving the aligned sequence concatenated to the reference sequence, mutation_type (insertion, deletion, or complex), and mutation length)
#' @export
#'
#' @examples
#' allele_table_fp <- "/Users/timbarry/research_offsite/external/bauer-lab/rhampseq_bcl11a/crispresso_output_1617/GUIDE0724-africa-edited/CRISPResso_on_1617_OT_0000/Alleles_frequency_table_around_sgRNA_CTAACAGTTGCTTTTATCAC.txt"
#' allele_feature_table <- process_crispresso_allele_table_cas9(allele_table_fp)
process_crispresso_allele_table_cas9 <- function(allele_table_fp) {
  allele_tab <- readr::read_delim(allele_table_fp, col_types = "ccliiidd")
  colnames(allele_tab) <- c("aligned_seq", "ref_seq", "unedited", "n_deleted", "n_inserted", "n_mutated", "n_reads", "percent_reads")

  # 1. separate mutated from unmutated alleles (as determined by CRISPResso)
  unmutated_tab <- allele_tab |> dplyr::filter(unedited)
  mutated_tab <- allele_tab |> dplyr::filter(!unedited)

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

  # 5. determine starting point of mutation (relative to cut site)
  aligned_seq_list <- strsplit(mutated_tab$aligned_seq, split = "")
  ref_seq_list <- strsplit(mutated_tab$ref_seq, split = "")
  window_size <- length(aligned_seq_list[[1]])
  quant_window_left_bdy <- window_size/2L - 1L
  quant_window_right_bdy <- window_size/2L + 2L
  cut_site <- window_size/2L
  leftmost_indel_pos <- sapply(X = seq(1L, nrow(mutated_tab)), FUN = function(i) {
    s <- if (mutation_types_new[i] == "insertion") {
      find_leftmost_indel_pos(s = ref_seq_list[[i]], quant_window_left_bdy = quant_window_left_bdy, quant_window_right_bdy = quant_window_right_bdy)
    } else if (mutation_types_new[i] == "deletion") {
      find_leftmost_indel_pos(s = aligned_seq_list[[i]], quant_window_left_bdy = quant_window_left_bdy, quant_window_right_bdy = quant_window_right_bdy)
    } else { # complex
      min(find_leftmost_indel_pos(s = aligned_seq_list[[i]], quant_window_left_bdy = quant_window_left_bdy, quant_window_right_bdy = quant_window_right_bdy),
          find_leftmost_indel_pos(s = ref_seq_list[[i]], quant_window_left_bdy = quant_window_left_bdy, quant_window_right_bdy = quant_window_right_bdy))
    }
    s - cut_site
  })

  # 6. prepare the output data frame
  out_df <- rbind(
    data.frame(aligned_seq = NA_character_, ref_seq = NA_character_,
               mutation_type = NA_character_,  mutation_length = NA,
               leftmost_indel_pos = NA, read_count = sum(unmutated_tab$n_reads),
               modified = FALSE),
    data.frame(aligned_seq = mutated_tab$aligned_seq, ref_seq = mutated_tab$ref_seq,
               mutation_type = mutation_types_new, mutation_length = mutation_lengths_new,
               leftmost_indel_pos = leftmost_indel_pos, read_count = mutated_tab$n_reads,
               modified = TRUE)
    )
  return(out_df)
}

find_leftmost_indel_pos <- function(s, quant_window_left_bdy, quant_window_right_bdy) {
  # find the contiguous hyphen substrings
  hyphen_idxs <- sort(grep(pattern = "-", fixed = TRUE, x = s))
  breaks <- c(TRUE, diff(hyphen_idxs) != 1)
  groups <- cumsum(breaks)
  n_groups <- max(groups)
  substring_idxs <- lapply(X = seq(1L, n_groups), FUN = function(i) hyphen_idxs[groups == i])

  # for each substring, determine (i) its leftmost position and (ii) whether it overlaps the quantification window
  substring_info <- sapply(X = substring_idxs, FUN = function(curr_substring) {
    overlaps_quant_window <- any(curr_substring >= quant_window_left_bdy & curr_substring <= quant_window_right_bdy)
    leftmost_pos <- min(curr_substring)
    c(overlaps_quant_window = overlaps_quant_window, leftmost_pos = leftmost_pos)
  }) |> t()
  leftmost_indel_pos <- min(substring_info[substring_info[,1L] == 1L, 2])
  return(leftmost_indel_pos)
}
