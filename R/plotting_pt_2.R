revlog_trans <- function(base = exp(1)) {
  trans <- function(x) {
    -log(x, base)
  }
  inv <- function(x) {
    base^(-x)
  }
  scales::trans_new(
    name = paste("revlog-", base, sep = ""),
    transform = trans,
    inverse = inv,
    breaks = scales::log_breaks(base = base),
    domain = c(1e-100, Inf)
  )
}

create_umi_hist_plot <- function(result_df, count_df, grna_spacer, pam, i) {
  ref_genome_fp <- paste0(.get_config_path("REF_GENOME_DIR"), "hg38_main_chroms.fa")
  fa <- Rsamtools::FaFile(ref_genome_fp)
  open(fa)
  x_range <- c(-30, 30)
  alignment_shift <- x_range[1] - 1
  curr_lead_base <- result_df[i,]$lead_base
  curr_chr <- as.character(seqnames(result_df[i,]))
  coord_range <- curr_lead_base + x_range
  count_df_shift <- count_df |>
    dplyr::filter(chr == curr_chr, coord >= coord_range[1], coord <= coord_range[2]) |>
    dplyr::mutate(coord = (coord - curr_lead_base))
  query_coord <- GenomicRanges::GRanges(curr_chr, IRanges(coord_range[1], coord_range[2]))
  query_seq <- Biostrings::getSeq(fa, query_coord)[[1]]
  grna_spacer_biostring <- Biostrings::DNAString(grna_spacer)
  pam_biostring <- Biostrings::DNAString(pam)
  alignment_res <- align_spacer_seq(query_seq = query_seq, grna_spacer = grna_spacer_biostring, pam_pattern = pam_biostring)
  if (nrow(alignment_res) >= 1L) {
    best_alignment_res <- alignment_res[1,]
    cut_base_range <- c(best_alignment_res[["cut_base_start"]], best_alignment_res[["cut_base_end"]]) + alignment_shift
  } else {
    cut_base_range <- grna_spacer_biostring <- NULL
    best_alignment_res <- list()
  }
  p <- make_local_scatterplot(count_df_sub = count_df_shift,
                              x_range = x_range,
                              log_trans = FALSE,
                              target_seq = query_seq,
                              grna_spacer = grna_spacer_biostring,
                              pam_strand = best_alignment_res[["pam_strand"]],
                              protospacer_range = c(best_alignment_res[["protospacer_start"]], best_alignment_res[["protospacer_end"]]) + alignment_shift,
                              pam_range = c(best_alignment_res[["pam_start"]], best_alignment_res[["pam_end"]]) + alignment_shift,
                              cut_base_range = cut_base_range)
    title <- paste0(paste0(seqnames(result_df[i,]), ":", ranges(result_df[i,])), "\n",
                  "gRNA alignment quality: ", best_alignment_res$score, "/20")
    p <- p + patchwork::plot_annotation(title)
  return(p)
}


#' Make local scatterplot
#'
#' @param count_df
#' @param x_range
#' @param log_trans
#' @param point_size
#' @param col
#' @param title
#' @param target_seq
#' @param grna_spacer
#' @param pam_strand
#' @param protospacer_range
#' @param pam_range
#' @param cut_base_range
#' @export
make_local_scatterplot <- function(count_df_sub, x_range = NULL, log_trans = FALSE, point_size = 1,
                                   col = c("dodgerblue3", "firebrick")[1], title = NULL, target_seq = NULL,
                                   grna_spacer = NULL, pam_strand = NULL, protospacer_range = NULL, pam_range = NULL,
                                   cut_base_range = NULL) {
  count_df_sub_plus <- count_df_sub |> dplyr::filter(strand == "+")
  count_df_sub_minus <- count_df_sub |> dplyr::filter(strand == "-")

  # construct label_df
  if (is(target_seq, "DNAString")) {
    target_seq <- strsplit(as.character(target_seq), split = "")[[1]]
  }
  if (is(grna_spacer, "DNAString")) {
    if (pam_strand == "-") grna_spacer <- Biostrings::reverseComplement(grna_spacer)
    grna_spacer <- strsplit(as.character(grna_spacer), split = "")[[1]]
  }
  label_df <- data.frame(x = seq(x_range[1], x_range[2])) |>
    dplyr::mutate(base = target_seq, y = -0.1) |>
    dplyr::mutate(base_type = "target")
  if (!is.null(grna_spacer)) {
    label_df <- label_df |>
      dplyr::mutate(pam_site = (x >= pam_range[1] & x <= pam_range[2])) |>
      dplyr::mutate(base_type = ifelse(pam_site, "pam", "target"),
                    pam_site = NULL)
    spacer_label_df <- data.frame(x = seq(protospacer_range[1], protospacer_range[2]),
                                  base = grna_spacer, y = 0.1, base_type = "spacer")
    label_df <- rbind(label_df, spacer_label_df)
  }

  # base plot
  make_base_plot <- function(curr_count_df_sub) {
    ggplot2::ggplot(data = curr_count_df_sub, mapping = ggplot2::aes(x = coord, y = count)) +
      ggplot2::geom_segment(ggplot2::aes(x = coord, xend = coord, y = if (log_trans) 1 else 0, yend = count)) +
      ggplot2::geom_point(size = point_size, col = col) +
      ggplot2::theme_bw(base_size = 10) + ggplot2::xlab("Coordinate") +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     axis.text.x  = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(0.0, 5.5, 0.0, 5.5)) +
      ggplot2::scale_x_continuous(limits = range(label_df$x))
  }
  p_plus <- make_base_plot(count_df_sub_plus)
  p_minus <- make_base_plot(count_df_sub_minus)

  # y-axis scale
  y_max <- max(c(count_df_sub_plus$count, count_df_sub_minus$count))
  y_limits <- if (log_trans) c(1, y_max) else c(0, y_max)
  if (log_trans) {
    p_plus <- p_plus +
      ggplot2::scale_y_continuous(trans = scales::log10_trans(),
                                  expand = c(0.01, 0),
                                  limits = y_limits) +
      ggplot2::ylab("")
    p_minus <- p_minus +
      ggplot2::scale_y_continuous(trans = revlog_trans(base = 10),
                                  expand = c(0.01, 0),
                                  limits = y_limits[2:1]) +
      ggplot2::ylab("")
  } else {
    p_plus <- p_plus + ggplot2::ylab("") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.01, add = 0),
                                  limits = y_limits)
    p_minus <- p_minus + ggplot2::ylab("") +
      ggplot2::scale_y_continuous(trans = scales::reverse_trans(),
                                  expand = ggplot2::expansion(mult = 0.01, add = 0),
                                  limits = y_limits[c(2L, 1L)])
  }

  p_middle <- ggplot2::ggplot() +
    ggplot2::geom_text(ggplot2::aes(x = x, y = y, label = base, col = base_type), data = label_df) +
    ggplot2::scale_color_manual(values = c("pam" = "red", "target" = "black", "spacer" = "darkorchid")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y  = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(0.0, 5.5, 0.0, 5.5),
                   legend.position = "none") +
    ggplot2::ylab(paste0("UMI count ", if (log_trans) "(log)" else "(linear)")) +
    ggplot2::scale_y_continuous(limits = c(-0.17, 0.15))

  if (!is.null(cut_base_range)) {
    cut_spot_x <- mean(cut_base_range)
    p_middle <- p_middle + ggplot2::geom_vline(xintercept = cut_spot_x, col = "orange")
  }

  library(patchwork)
  p_all <- (p_plus /  p_middle / p_minus) +
    plot_layout(heights = c(1, 0.13, 1), axes = "collect")

  return(p_all)
}


make_histogram_plot <- function(y, fit, result_df, model) {
  # compute fitted density
  x_range <- seq(0, max(y))
  if (model == "nb") {
    expected <- dnbinom(x = x_range, size = fit[["theta"]], mu = fit[["mu"]]) * length(y)
  } else if (model == "poisson") {
    expected <- dpois(x = x_range, lambda = fit[["mu"]]) * length(y)
  } else {
    stop("Model not recognized.")
  }
  shifted_density_df <- data.frame(count = x_range, expected = expected)

  # shift
  x_range <- x_range + 1
  y <- y + 1
  shifted_density_df$count <- shifted_density_df$count + 1

  # create untrans histogram
  p_model_untrans <- ggplot2::ggplot(data = data.frame(y = y), mapping = ggplot2::aes(x = y)) +
    ggplot2::geom_histogram(binwidth = 1, col = "black", fill = "white") +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
      data = shifted_density_df,
      ggplot2::aes(x = count, y = expected),
      linewidth = 0.9,
      color = "firebrick") +
    ggplot2::ylab("Frequency") +
    ggplot2::ggtitle("Linear y-axis") + ggplot2::xlab("UMI count")

  # check if discoveries are present; if so, draw line
  if (any(result_df$significant_hit)) {
    rejection_thresh <- (result_df |> plyranges::filter(significant_hit))$umi_count |> min()
    p_model_untrans <- p_model_untrans + ggplot2::geom_vline(xintercept = rejection_thresh, col = "blue", linetype = "dashed")
  }
  p_model_trans <- p_model_untrans +
    ggplot2::scale_y_continuous(transform = scales::pseudo_log_trans(), breaks = c(0, 10^seq(0, 5))) +
    ggplot2::ggtitle("Log y-axis")
  out <- list(plot_untrans = p_model_untrans, plot_trans = p_model_trans)
  return(out)
}


#' Plot read processing results
#'
#' @param read_proc_df the output of `run_read_qc_on_sample()`
#'
#' @returns
#' @export
#'
#' @examples
#' sample_dir <- "/Users/timbarry/research_code/genethoff-nf/demo/results/Jing_AAVS1_n1-10_Donor-Seq_AAVS1_GSP_plus_1"
#' read_proc_df <- run_read_qc_on_sample(sample_dir)
#' plot_read_processing_results(read_proc_df)
plot_read_processing_results <- function(read_proc_df) {
  # paired-end plot
  to_plot_paired <- read_proc_df |>
    dplyr::filter(category %in% c("n_reads_missing_dsodn_tag", "n_reads_too_short",
                                  "n_reads_unaligned", "n_reads_poorly_aligned",
                                  "n_reads_multimapped", "n_reads_good_alignment")) |>
    dplyr::mutate(Category = factor(x = category,
                                    levels = c("n_reads_missing_dsodn_tag", "n_reads_too_short", "n_reads_unaligned", "n_reads_poorly_aligned", "n_reads_multimapped", "n_reads_good_alignment") ,
                                    labels = c("N reads missing dsODN tag", "N reads too short", "N reads unaligned", "N reads poorly aligned", "N reads multimapped", "N reads w/ good alignment")))
  p_paired <- ggplot2::ggplot(to_plot_paired, ggplot2::aes(x = stage, y = n_reads, fill = Category)) +
    ggplot2::geom_bar(stat = "identity") + ggplot2::theme_bw()

  # r2 plot
  to_plot_r2 <- read_proc_df |>
    dplyr::filter(category %in% c("n_reads_too_short_r2", "n_reads_unaligned_r2",
                                  "n_reads_poorly_aligned_r2", "n_reads_multimapped_r2", "n_reads_good_alignment_r2")) |>
    dplyr::mutate(Category = factor(x = category,
                                    levels = c("n_reads_too_short_r2", "n_reads_unaligned_r2", "n_reads_poorly_aligned_r2", "n_reads_multimapped_r2", "n_reads_good_alignment_r2") ,
                                    labels = c("N reads too short (R2)", "N reads unaligned (R2)", "N reads poorly aligned (R2)", "N reads multimapped (R2)", "N reads w/ good alignment (R2)")))

  p_r2 <- ggplot2::ggplot(data = to_plot_r2, mapping = ggplot2::aes(x = stage, y = n_reads, fill = Category)) +
    ggplot2::geom_bar(stat = "identity") + ggplot2::theme_bw()

  # combined plot
  to_plot_combined <- read_proc_df |>
    dplyr::filter(category %in% c("n_reads_missing_dsodn_tag", "n_reads_poorly_aligned", "n_reads_multimapped", "n_reads_good_alignment",
                                  "n_reads_too_short_r2", "n_reads_unaligned_r2", "n_reads_poorly_aligned_r2", "n_reads_multimapped_r2", "n_reads_good_alignment_r2")) |>
    dplyr::mutate(Category = factor(x = category,
                                    levels = c("n_reads_missing_dsodn_tag", "n_reads_poorly_aligned", "n_reads_multimapped", "n_reads_good_alignment", "n_reads_too_short_r2", "n_reads_unaligned_r2", "n_reads_poorly_aligned_r2", "n_reads_multimapped_r2", "n_reads_good_alignment_r2") ,
                                    labels = c("N reads missing dsODN tag", "N reads poorly aligned (paired-end)", "N reads multimapped (paired-end)", "N reads w/ good alignment (paired-end)", "N reads too short (R2)", "N reads unaligned (R2)", "N reads poorly aligned (R2)", "N reads multimapped (R2)", "N reads w/ good alignment (R2)")))

  p_combined <- ggplot2::ggplot(data = to_plot_combined, mapping = ggplot2::aes(x = "Sample", y = n_reads, fill = Category)) +
    ggplot2::geom_bar(stat = "identity") + ggplot2::theme_bw() + ggplot2::ylab("N reads") + ggplot2::xlab("")
}
