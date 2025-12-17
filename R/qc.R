fetch_number <- function(sample_dir, file_name, n_skip, posit_from_end = NA, posit_from_start = NA) {
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


get_alignment_qc <- function(sample_dir, align_log, align_qc_log) {
  n_reads_0_alignments <- fetch_number(sample_dir = sample_dir, file_name = align_log, n_skip = 1L, posit_from_start = 2L)
  n_reads_1_alignment <- fetch_number(sample_dir = sample_dir, file_name = align_log, n_skip = 2L, posit_from_start = 2L)
  n_reads_2_plus_alignments <- fetch_number(sample_dir = sample_dir, file_name = align_log, n_skip = 3L, posit_from_start = 2L)

  n_reads_good_mapq <- fetch_number(sample_dir = sample_dir, file_name = align_qc_log, n_skip = 0L, posit_from_end = 0L)
  n_reads_1_align_good_mapq <- fetch_number(sample_dir = sample_dir, file_name = align_qc_log, n_skip = 1L, posit_from_end = 0L)
  n_reads_1_alignment_poor_mapq <- n_reads_1_alignment - n_reads_1_align_good_mapq
  n_reads_2_plus_align_good_mapq <- fetch_number(sample_dir = sample_dir, file_name = align_qc_log, n_skip = 2L, posit_from_end = 0L)
  n_reads_2_plus_align_poor_mapq <- n_reads_2_plus_alignments - n_reads_2_plus_align_good_mapq

  n_reads_unmapped_or_poorly_mapped <- n_reads_0_alignments + n_reads_1_alignment_poor_mapq
  n_reads_multimapped <- n_reads_2_plus_align_poor_mapq

  c(n_reads_good_mapq = n_reads_good_mapq,
    n_reads_0_alignments = n_reads_0_alignments,
    n_reads_unmapped_or_poorly_mapped = n_reads_unmapped_or_poorly_mapped,
    n_reads_poorly_mapped = n_reads_unmapped_or_poorly_mapped - n_reads_0_alignments,
    n_reads_multimapped = n_reads_multimapped)
}


#' Run read QC on sample
#'
#' @param sample_dir
#'
#' @returns
#' @export
#'
#' @examples
#' sample_dirs <- paste0("/Users/timbarry/research_code/genethoff-nf/demo/results/",
#' c("Jing_AAVS1_n1-10_Donor-Seq_AAVS1_GSP_plus_1", "Jing_AAVS1_n1-10_Donor-Seq_IL2RG_GSP_minus_1"))
#' qc_df <- run_read_qc_on_sample(sample_dirs)
#'
run_read_qc_on_sample <- function(sample_dirs) {
  df_out <- lapply(X = sample_dirs, FUN = function(sample_dir) {
    # 0. total number of reads
    n_total_reads <- fetch_number(sample_dir = sample_dir, file_name = "/1_trim_5_prime_tag.log", n_skip = 3L, posit_from_end = 0L)

    # 1. number of reads filtered out due to missing dsodn tag
    n_reads_missing_tag <- fetch_number(sample_dir = sample_dir, file_name = "/1_trim_5_prime_tag.log", n_skip = 8L, posit_from_end = 1L)

    # 3. paired-end reads filtered out as too short
    n_reads_too_short_paired_end <- fetch_number(sample_dir = sample_dir, file_name = "/3_filter_by_length.log", n_skip = 7L, posit_from_end = 1L)

    # 4. get pairwise alignment metrics
    paired_align_metrics <- get_alignment_qc(sample_dir = sample_dir, align_log = "/4_align_paired_end_reads.log", align_qc_log = "/5_paired_end_alignment_qc.log")

    # 5. number of r2 reads at start of rescue step
    n_r2_reads_rescue_start <- paired_align_metrics[["n_reads_0_alignments"]] + n_reads_too_short_paired_end
    n_reads_too_short_r2 <- fetch_number(sample_dir = sample_dir, file_name = "/6_filter_r2_leftovers_by_length.log", n_skip = 7L, posit_from_end = 1L)

    # 6. get r2 alignment metrics
    r2_alignment_metrics <- get_alignment_qc(sample_dir = sample_dir, align_log = "/7_align_r2_reads.log", align_qc_log = "/8_r2_alignment_qc.log")

    # construct output
    ret <- data.frame("n_total_reads" = n_total_reads,
                      "n_reads_missing_dsodn_tag" = n_reads_missing_tag,
                      "n_reads_too_short" = n_reads_too_short_paired_end,
                      "n_reads_unaligned" = paired_align_metrics[["n_reads_0_alignments"]],
                      "n_reads_poorly_aligned" = paired_align_metrics[["n_reads_poorly_mapped"]],
                      "n_reads_multimapped" = paired_align_metrics[["n_reads_multimapped"]],
                      "n_reads_good_alignment" = paired_align_metrics[["n_reads_good_mapq"]],
                      "n_reads_r2_rescue_start" = paired_align_metrics[["n_reads_0_alignments"]] + n_reads_too_short_paired_end,
                      "n_reads_too_short_r2" = n_reads_too_short_r2,
                      "n_reads_unaligned_r2" = r2_alignment_metrics[["n_reads_0_alignments"]],
                      "n_reads_poorly_aligned_r2" = r2_alignment_metrics[["n_reads_poorly_mapped"]],
                      "n_reads_multimapped_r2" = r2_alignment_metrics[["n_reads_multimapped"]],
                      "n_reads_good_alignment_r2" = r2_alignment_metrics[["n_reads_good_mapq"]])
      tidyr::pivot_longer(data = ret, cols = colnames(ret)) |>
      setNames(c("category", "n_reads")) |>
      dplyr::mutate(sample = BiocGenerics::basename(sample_dir))
  }) |> data.table::rbindlist() |> dplyr::mutate(sample = factor(sample))
  return(df_out)
}
