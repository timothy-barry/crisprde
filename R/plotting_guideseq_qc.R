#' Plot read processing results
#'
#' @param qc_df the output of `run_read_qc_on_sample()`
#' @param plot_type type of plot, either "paired_end", "r2", or "combined"
#' @param normalize_counts logical indicating whether to normalize each sample's
#'   displayed categories to sum to one.
#'
#' @returns a `ggplot2` stacked bar plot.
#'
#' @details
#' The read-processing plots use the disjoint QC categories returned by
#' `run_read_qc_on_sample()`. Multimappers are shown as separate discarded and
#' retained categories. The retained multimapper category corresponds to
#' exact-tie multimappers that passed the pipeline's post-alignment QC when
#' multimapper retention was enabled.
#' @export
#'
#' @examples
#' sample_dirs <- paste0("/Users/timbarry/research_code/genethoff-nf/demo/results/",
#' c("Jing_AAVS1_n1-10_Donor-Seq_AAVS1_GSP_plus_1", "Jing_AAVS1_n1-10_Donor-Seq_IL2RG_GSP_minus_1"))
#' qc_df <- run_read_qc_on_sample(sample_dirs)
#' plot_read_processing_results(qc_df = qc_df, plot_type = "combined")
plot_read_processing_results <- function(qc_df, plot_type = "combined", normalize_counts = FALSE) {
  if (plot_type == "paired_end") {
    # paired-end plot
    plot_levels <- c("n_reads_missing_dsodn_tag", "n_reads_too_short", "n_reads_unaligned", "n_reads_poorly_aligned", "n_reads_discarded_multimapped", "n_reads_retained_multimapped", "n_reads_good_alignment")
    plot_labels <- c("N reads missing dsODN tag", "N reads too short", "N reads unaligned", "N reads poorly aligned", "N reads discarded multimapped", "N reads retained multimapped", "N reads w/ good alignment")
    plot_colors <- c("#7FFFD4", "#8B2323", "#87CEFA", "#7A378B", "#F0E68C", "#FFB000", "#008B00")
  } else if (plot_type == "r2") {
    # r2 plot
    plot_levels <- c("n_reads_too_short_r2", "n_reads_unaligned_r2", "n_reads_poorly_aligned_r2", "n_reads_discarded_multimapped_r2", "n_reads_retained_multimapped_r2", "n_reads_good_alignment_r2")
    plot_labels <- c("N reads too short (R2)", "N reads unaligned (R2)", "N reads poorly aligned (R2)", "N reads discarded multimapped (R2)", "N reads retained multimapped (R2)", "N reads w/ good alignment (R2)")
    plot_colors <- c("#8B2323","#87CEFA", "#7A378B", "#F0E68C", "#FFB000", "#008B00")
  } else if (plot_type == "combined") {
    # combined plot
    plot_levels <- c("n_reads_missing_dsodn_tag", "n_reads_poorly_aligned", "n_reads_discarded_multimapped", "n_reads_retained_multimapped", "n_reads_good_alignment", "n_reads_too_short_r2", "n_reads_unaligned_r2", "n_reads_poorly_aligned_r2", "n_reads_discarded_multimapped_r2", "n_reads_retained_multimapped_r2", "n_reads_good_alignment_r2")
    plot_labels <- c("N reads missing dsODN tag", "N reads poorly aligned (paired-end)", "N reads discarded multimapped (paired-end)", "N reads retained multimapped (paired-end)", "N reads w/ good alignment (paired-end)", "N reads too short (R2)", "N reads unaligned (R2)", "N reads poorly aligned (R2)", "N reads discarded multimapped (R2)", "N reads retained multimapped (R2)", "N reads w/ good alignment (R2)")
    plot_colors <- c("#7FFFD4", "#7A378B", "#F0E68C", "#FFB000", "#008B00", "#8B2323", "#87CEFA", "#8B5742", "#EEE9E9", "#E69F00", "#104E8B")
  }
  to_plot <- qc_df |>
    dplyr::filter(category %in% plot_levels) |>
    dplyr::mutate(Category = factor(x = category, levels = plot_levels, labels = plot_labels))
  if (normalize_counts) {
    to_plot <- to_plot |> dplyr::group_by(sample) |> dplyr::mutate(n_reads = n_reads/sum(n_reads))
  }

  p <- ggplot2::ggplot(to_plot, ggplot2::aes(x = sample, y = n_reads, fill = Category)) +
    ggplot2::geom_bar(stat = "identity") + ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = plot_colors) +
    ggplot2::ylab("N reads") + ggplot2::xlab("Sample") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  return(p)
}
