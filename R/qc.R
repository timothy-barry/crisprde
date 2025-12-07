#' Run read QC on sample
#'
#' @param sample_dir
#'
#' @returns
#' @export
#'
#' @examples
#' sample_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"),
#' "guideseq_elane/count_tables_no_multimap/CD34-WT-Cas9-ELANE-e3SA-GSPneg-5uM-S89_S97/")
run_read_qc_on_sample <- function(sample_dir) {
  fetch_number <- function(file_name, n_skip, posit_from_end = NA, posit_from_start = NA) {
    df <- read.delim(paste0(sample_dir, file_name),
                            nrows = 1, skip = n_skip, header = TRUE, col.names = "n_reads")
    if (is.integer(df$n_reads)) {
      out <- df$n_reads
    } else {
      str_split_out <- strsplit(x = df$n_reads, split = "\\s+")[[1]]
      if (!is.na(posit_from_end)) {
        n_reads <- str_split_out[length(str_split_out) - posit_from_end]
      } else {
        n_reads <- str_split_out[posit_from_start]
      }
      out <- gsub(pattern = ",", replacement = "", x = n_reads) |> as.integer()
    }
    return(out)
  }

  # 0. total number of reads
  n_total_reads <- fetch_number(file_name = "/1_trim_5_prime_tag.log", n_skip = 3L, posit_from_end = 0L)

  # 1. number of reads filtered out due to missing dsodn tag
  n_reads_missing_tag <- fetch_number(file_name = "/1_trim_5_prime_tag.log", n_skip = 8L, posit_from_end = 1L)

  # 3. paired-end reads filtered out as too short
  n_reads_too_short_paired_end <- fetch_number(file_name = "/3_filter_by_length.log", n_skip = 7L, posit_from_end = 1L)

  # 4. align paired end reads qc
  n_reads_0_alignments <- fetch_number(file_name = "/4_align_paired_end_reads.log", n_skip = 1L, posit_from_start = 2L)
  n_reads_1_alignment <- fetch_number(file_name = "/4_align_paired_end_reads.log", n_skip = 2L, posit_from_start = 2L)
  n_reads_1_plus_alignments <- fetch_number(file_name = "/4_align_paired_end_reads.log", n_skip = 3L, posit_from_start = 2L)
  n_reads_passing_mapq <- fetch_number(file_name = "/5_paired_end_alignment_qc.log", n_skip = 0L, posit_from_start = 1L)
  n_reads_removed_multimapping <- n_reads_1_plus_alignments - (n_reads_passing_mapq - n_reads_1_alignment)

  # 5. n r2 reads at start of rescue step
  n_r2_reads_rescue_start <- n_reads_0_alignments + n_reads_too_short_paired_end
  n_reads_too_short_r2 <- fetch_number(file_name = "/6_filter_r2_leftovers_by_length.log", n_skip = 7L, posit_from_end = 1L)

  # 6. align r2 reads
  n_reads_0_alignments_r2 <- fetch_number(file_name = "/7_align_r2_reads.log", n_skip = 1L, posit_from_start = 2L)
  n_reads_1_alignment_r2 <- fetch_number(file_name = "/7_align_r2_reads.log", n_skip = 2L, posit_from_start = 2L)
  n_reads_1_plus_alignments_r2 <- fetch_number(file_name = "/7_align_r2_reads.log", n_skip = 3L, posit_from_start = 2L)
  n_reads_passing_mapq_r2 <- fetch_number(file_name = "/8_r2_alignment_qc.log", n_skip = 0L, posit_from_start = 1L)
  n_reads_removed_multimapping_r2 <- n_reads_1_plus_alignments_r2 - (n_reads_passing_mapq_r2 - n_reads_1_alignment_r2)

  # construct output
  ret <- data.frame("n_total_reads" = n_total_reads,
           "n_reads_missing_dsodn_tag" = n_reads_missing_tag,
           "n_reads_too_short" = n_reads_too_short_paired_end,
           "n_reads_unaligned" = n_reads_0_alignments,
           "n_reads_poor_mapq_or_multimap" = n_reads_removed_multimapping,
           "n_reads_good_alignment" = n_reads_passing_mapq,
           "n_reads_r2_rescue_start" = n_r2_reads_rescue_start,
           "n_reads_too_short_r2" = n_reads_too_short_r2,
           "n_reads_unaligned_r2" = n_reads_0_alignments_r2,
           "n_reads_poor_mapq_or_multimap_r2" = n_reads_removed_multimapping_r2,
           "n_reads_good_alignment_r2" = n_reads_passing_mapq_r2)
  tidyr::pivot_longer(ret, cols = colnames(ret)) |>
    setNames(c("category", "n_reads")) |>
    dplyr::mutate(stage = ifelse(grepl(pattern = "r2", x = category), "r2", "paired_end"))
}
