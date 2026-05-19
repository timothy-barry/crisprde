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

  n_reads_passing_qc <- fetch_number(sample_dir = sample_dir, file_name = align_qc_log, n_skip = 0L, posit_from_end = 0L)
  n_reads_1_align_good_mapq <- fetch_number(sample_dir = sample_dir, file_name = align_qc_log, n_skip = 1L, posit_from_end = 0L)
  n_reads_1_alignment_poor_mapq <- n_reads_1_alignment - n_reads_1_align_good_mapq
  n_reads_2_plus_align_good_mapq <- fetch_number(sample_dir = sample_dir, file_name = align_qc_log, n_skip = 2L, posit_from_end = 0L)

  n_reads_good_alignment <- n_reads_1_align_good_mapq + n_reads_2_plus_align_good_mapq
  n_reads_retained_multimapped <- n_reads_passing_qc - n_reads_good_alignment
  n_reads_discarded_multimapped <- n_reads_2_plus_alignments -
    n_reads_2_plus_align_good_mapq -
    n_reads_retained_multimapped

  counts_to_check <- c(n_reads_1_alignment_poor_mapq,
                       n_reads_retained_multimapped,
                       n_reads_discarded_multimapped)
  if (any(is.na(counts_to_check)) || any(counts_to_check < 0L)) {
    stop("Alignment QC counts are inconsistent for sample directory: ", sample_dir,
         call. = FALSE)
  }

  c(n_reads_good_alignment = n_reads_good_alignment,
    n_reads_0_alignments = n_reads_0_alignments,
    n_reads_poorly_mapped = n_reads_1_alignment_poor_mapq,
    n_reads_discarded_multimapped = n_reads_discarded_multimapped,
    n_reads_retained_multimapped = n_reads_retained_multimapped)
}


#' Run read QC on sample
#'
#' @param sample_dirs character vector of GENETHOFF/Donor-seq sample output
#'   directories.
#'
#' @returns a data frame of read counts by QC category and sample. Multimapping
#'   reads are split into discarded multimappers and retained exact-tie
#'   multimappers.
#'
#' @details
#' The output has one row per `(sample, category)` pair with columns `category`,
#' `n_reads`, and `sample`. For each alignment step, `n_reads_good_alignment`
#' counts reads passing the high-MAPQ criteria, while
#' `n_reads_retained_multimapped` counts low-MAPQ exact-tie multimappers retained
#' by the GENETHOFF/Donor-seq pipeline when multimapper retention is enabled.
#' `n_reads_discarded_multimapped` counts remaining multimappers that did not
#' pass post-alignment QC. The corresponding R2 rescue categories use the `_r2`
#' suffix. The function infers these buckets from the pipeline logs; callers do
#' not need to specify whether `keep_multimapped_reads` was `TRUE` or `FALSE`.
#'
#' In a stacked read-processing plot, the paired-end categories
#' `n_reads_missing_dsodn_tag`, `n_reads_too_short`, `n_reads_unaligned`,
#' `n_reads_poorly_aligned`, `n_reads_discarded_multimapped`,
#' `n_reads_retained_multimapped`, and `n_reads_good_alignment` should sum to
#' `n_total_reads`. The analogous R2 rescue categories should sum to
#' `n_reads_r2_rescue_start`.
#' @export
#'
#' @examples
#' sample_dirs <- paste0("/Users/timbarry/research_offsite/external/bauer-lab/guideseq_sbds/count_tables_with_multimap/",
#' c("Jing_Max112625_1_SBDSP1_minus", "Jing_Max112625_1_SBDSP1_plus"))
#' out <- run_read_qc_on_sample(sample_dirs)
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
                      "n_reads_discarded_multimapped" = paired_align_metrics[["n_reads_discarded_multimapped"]],
                      "n_reads_retained_multimapped" = paired_align_metrics[["n_reads_retained_multimapped"]],
                      "n_reads_good_alignment" = paired_align_metrics[["n_reads_good_alignment"]],
                      "n_reads_r2_rescue_start" = n_r2_reads_rescue_start,
                      "n_reads_too_short_r2" = n_reads_too_short_r2,
                      "n_reads_unaligned_r2" = r2_alignment_metrics[["n_reads_0_alignments"]],
                      "n_reads_poorly_aligned_r2" = r2_alignment_metrics[["n_reads_poorly_mapped"]],
                      "n_reads_discarded_multimapped_r2" = r2_alignment_metrics[["n_reads_discarded_multimapped"]],
                      "n_reads_retained_multimapped_r2" = r2_alignment_metrics[["n_reads_retained_multimapped"]],
                      "n_reads_good_alignment_r2" = r2_alignment_metrics[["n_reads_good_alignment"]])
      tidyr::pivot_longer(data = ret, cols = colnames(ret)) |>
      setNames(c("category", "n_reads")) |>
      dplyr::mutate(sample = BiocGenerics::basename(sample_dir))
  }) |> data.table::rbindlist() |> dplyr::mutate(sample = factor(sample))
  return(df_out)
}
