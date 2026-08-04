#' Make GUIDE-seq QQ plot
#'
#' @param res_df a data frame with column `p_value` giving the p-value of a given locus
#'
#' @returns a qq plot of the results
#' @export
#'
#' @examples
#' set.seed(7)
#' # NULL DATA
#' pi <- c(0.05, 0.1, 0.02)
#' mu_vect <- c(10, 6, 15)
#' theta_vect <- c(2, 5, 0.3)
#' m <- 10000
#' null_dat <- simulate_multirep_guideseq_data(pi, mu_vect, theta_vect, m)
#'
#' # ALTERNATIVE DATA
#' pi <- c(0.5, 0.8, 0.6)
#' mu_vect <- c(250, 200, 50)
#' theta_vect <- c(50, 50, 60)
#' m_alt <- 15
#' alt_dat <- simulate_multirep_guideseq_data(pi, mu_vect, theta_vect, m_alt)
#'
#' # COMBINED DATA
#' Y_mat <- cbind(alt_dat, null_dat)
#' colnames(Y_mat) <- paste0("window_", seq_len(ncol(Y_mat)))
#' incorporate_occupancy_info <- TRUE
#' multiplicity_alpha <- 0.2
#'
#' # RUN METHOD
#' res_df <- run_multireplicate_guideseq_method(Y_mat, incorporate_occupancy_info = TRUE, multiplicity_alpha = 0.2)
#' res_df$true_editing <- c(rep(TRUE, m_alt), rep(FALSE, m))
#' p_all <- make_guideseq_qq_plot(res_df, color_ground_truth = TRUE)
#' p_null <- make_guideseq_qq_plot(res_df[seq(m_alt+1, nrow(res_df)),], color_ground_truth = FALSE)
make_guideseq_qq_plot <- function(res_df, color_ground_truth = FALSE, rev_log_trans = TRUE,
                                  min_p = 1e-16, annotate_discoveries = FALSE, point_col = "black",
                                  rej_threshold_col = "blue", annotation_size = 1.5, annotation_color = "red") {
  # restrict to minimum p-value, get rejection threshold
  if (!is.na(min_p)) {
    res_df <- res_df |> dplyr::mutate(p_value = ifelse(p_value < min_p, min_p, p_value))
  }
  rejected_p_vals <- res_df |>
    dplyr::filter(nominated_window) |>
    dplyr::pull(p_value)
  rejection_threshold <- if (length(rejected_p_vals) == 0L) NULL else max(rejected_p_vals)

  if (color_ground_truth) {
    res_df <- res_df |> dplyr::mutate(true_editing = ifelse(true_editing, "Truly edited", "Truly unedited"))
    mapping <- ggplot2::aes(y = p_value, col = true_editing, size = true_editing, group = 1L)
  } else {
    mapping <- ggplot2::aes(y = p_value)
  }
  p <- ggplot2::ggplot(data = res_df, mapping = mapping) +
    stat_qq_band() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black")
  if (rev_log_trans) {
    p <- p + ggplot2::scale_x_continuous(trans = revlog_trans(base = 10)) +
      ggplot2::scale_y_continuous(trans = revlog_trans(base = 10), limits = c(1, min_p))
  } else {
   p <- p + ggplot2::scale_x_continuous(trans = scales::reverse_trans()) +
     ggplot2::scale_y_continuous(trans = scales::reverse_trans(), limits = c(1, min_p))
  }
  if (color_ground_truth) {
    p <- p + stat_qq_points(ymin = min_p) +
      ggplot2::scale_color_manual(values = c("firebrick1", "black")) +
      ggplot2::labs(color = "Ground truth") +
      ggplot2::scale_size_manual(values = c("Truly edited" = 1.5, "Truly unedited" = 0.5)) +
      ggplot2::guides(size = "none")
  } else {
    p <- p + stat_qq_points(size = 0.8, ymin = min_p, col = point_col)
  }
  p <- p + ggplot2::geom_hline(yintercept = rejection_threshold, col = rej_threshold_col, linetype = "dashed")

  # annotate
  if (annotate_discoveries) {
    label_df <- res_df |> dplyr::mutate(
      p_value_plot = pmax(p_value, min_p),
      qq_x = stats::qunif(stats::ppoints(dplyr::n())[rank(p_value_plot, ties.method = "first")]),
      label_y = pmin(1, p_value_plot * 1.5)
    ) |>
      dplyr::filter(nominated_window)
    p <- p + ggplot2::geom_text(data = label_df,
                                mapping = ggplot2::aes(x = qq_x, y = label_y, label = window_label),
                                inherit.aes = FALSE, angle = 90,
                                color = annotation_color, size = annotation_size, hjust = 1, vjust = 0.5)
  }

  return(p)
}



#' Make local scatterplot
#'
#' @param annotated_df_sub the annotated data frame for one window. gRNA and DNA sequence information must be present.
#' @param title plot title
#'
#' @returns a ggplot of the local UMI count distribution
#' @export
#'
#' @examples
#' caliper_res <- readRDS("/Users/timbarry/research_offsite/projects/crisprde-project/guideseq/hyperparam_res_list.rds")
#'
#'
#' # plus strand PAM
#' res_df <- caliper_res$elane_cd34_wtcas9_e3sa$tuning_res$selected_trt_run$res_df
#' window_id <- res_df |> dplyr::slice(1L) |> dplyr::pull(window)
#' annotated_df_sub <- caliper_res$elane_cd34_wtcas9_e3sa$annotated_clustered_count_df_trt |> dplyr::filter(window == window_id)
#' p <- make_local_scatterplot(annotated_df_sub)
#'
#' # minus strand PAM
#' res_df <- caliper_res$bcl11a_293t_1620_sprycas9$tuning_res$selected_trt_run$res_df
#' window_id <- res_df |> dplyr::slice(1L) |> dplyr::pull(window)
#' annotated_df_sub <- caliper_res$bcl11a_293t_1620_sprycas9$annotated_clustered_count_df_trt |> dplyr::filter(window == window_id)
#' p <- make_local_scatterplot(annotated_df_sub)
make_local_scatterplot <- function(annotated_df_sub, title = NULL) {
  library(patchwork)

  count_df_sub_plus <- annotated_df_sub |> dplyr::filter(strand == "+")
  count_df_sub_minus <- annotated_df_sub |> dplyr::filter(strand == "-")

  # prepare sequences
  dna_seq <- (annotated_df_sub$homology_dna[1] |> strsplit(split = ""))[[1]]
  grna_spacer <- (annotated_df_sub$homology_gRNA[1] |> strsplit(split = ""))[[1]]
  pam_site <- seq(length(grna_spacer) - 2L, length(grna_spacer))
  pam_strand <- annotated_df_sub$homology_strand[1]
  grna_spacer[pam_site] <- ""
  grna_spacer[grna_spacer == "T"] <- "U"
  homology_start <- annotated_df_sub$homology_posit[1]
  if (pam_strand == "+") {
    x_range <- seq(homology_start, homology_start + length(dna_seq) - 1L) + 1L
  } else {
    x_range <- seq(homology_start + length(dna_seq), homology_start + 1L)
  }
  label_df <- data.frame(coord = x_range,
                         dna_seq = dna_seq,
                         grna_spacer = grna_spacer)
  label_df$base_type <- "protospacer"
  label_df$base_type[pam_site] <- "pam"
  cut_start_posit <- annotated_df_sub$homology_cut_start[1]
  cut_end_posit <- annotated_df_sub$homology_cut_end[1]

  # sum across replicates and primer orientations
  count_df_sub_aggr <- annotated_df_sub |>
    dplyr::select(coord, umi_count, strand) |>
    dplyr::group_by(coord, strand) |>
    dplyr::summarize(umi_count = sum(umi_count)) |> dplyr::ungroup()

  # create the df to plot
  create_df_to_plot <- function(label_df, count_df_sub_aggr, curr_strand) {
    df_to_plot <- dplyr::left_join(label_df,
                                   count_df_sub_aggr |>
                                     dplyr::filter(strand == curr_strand) |>
                                     dplyr::select(-strand),
                                   by = "coord")
  }
  minus_df_to_plot <- create_df_to_plot(label_df = label_df, count_df_sub_aggr = count_df_sub_aggr, curr_strand = "-")
  plus_df_to_plot <- create_df_to_plot(label_df = label_df, count_df_sub_aggr = count_df_sub_aggr, curr_strand = "+")

  make_base_plot <- function(curr_count_df_sub) {
    p <- ggplot2::ggplot(data = curr_count_df_sub |> na.omit(),
                         mapping = ggplot2::aes(x = coord, y = umi_count)) +
      ggplot2::geom_segment(ggplot2::aes(x = coord, xend = coord, y = 1, yend = umi_count)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw(base_size = 10) + ggplot2::xlab("Coordinate") +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     axis.text.x  = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(0.0, 5.5, 0.0, 5.5)) +
      ggplot2::scale_x_continuous(limits = range(label_df$coord))
  }
  p_plus <- make_base_plot(plus_df_to_plot)
  p_minus <- make_base_plot(minus_df_to_plot)

  # make top and bottom umi plots
  y_max <- max(c(plus_df_to_plot$umi_count, minus_df_to_plot$umi_count), na.rm = TRUE)
  y_limits <- c(0, y_max)
  p_plus <- p_plus + ggplot2::ylab("") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.025, add = 0),
                                limits = y_limits)
  p_minus <- p_minus + ggplot2::ylab("") +
    ggplot2::scale_y_continuous(trans = scales::reverse_trans(),
                                expand = ggplot2::expansion(mult = 0.025, add = 0),
                                limits = y_limits[c(2L, 1L)])

  # make middle plot
  to_plot <- tidyr::pivot_longer(data = label_df, cols = c("dna_seq", "grna_spacer"),
                      names_to = "dna_or_rna", values_to = "base_value") |>
    dplyr::mutate(y = ifelse(dna_or_rna == "dna_seq", 0, 0.04),
                  base_type_dna_or_rna = paste0(base_type, "_", dna_or_rna))

  p_middle <- ggplot2::ggplot() +
    ggplot2::geom_text(ggplot2::aes(x = coord, y = y, label = base_value, col = base_type_dna_or_rna),
                       data = to_plot, size = 3.2) +
    ggplot2::scale_color_manual(values = c("pam_dna_seq" = "red",
                                           "protospacer_dna_seq" = "black",
                                           "protospacer_grna_spacer" = "dodgerblue3")) +
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
      ggplot2::ylab("") +
      ggplot2::scale_y_continuous(limits = c(-0.02, 0.06),
                                  expand = ggplot2::expansion(mult = 0, add = 0))
  if (is.null(title)) {
    title <- paste0(annotated_df_sub$window[1], " (", pam_strand, " strand PAM)")
  }

      cut_spot_x <- mean(c(cut_start_posit, cut_end_posit))
      p_middle <- p_middle + ggplot2::geom_vline(xintercept = cut_spot_x, col = "orange")
  p_all <- (p_plus /  p_middle / p_minus) +
    plot_layout(heights = c(1, 0.22, 1), axes = "collect") +
    patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 10, hjust = 0.5))
    )

  return(p_all)
}
