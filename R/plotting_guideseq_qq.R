#' Make GUIDE-seq QQ plot
#'
#' @param res_df a data frame with column `p_value` giving the p-value of a given locus
#'
#' @returns a qq plot of the results
#' @export
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



#' Make a simple GUIDE-seq QQ plot
#'
#' @param p_values Numeric vector of p-values to compare with a uniform null
#'   distribution.
#' @param rev_log_trans Use negative log10 scales if TRUE (the default), or
#'   reversed linear scales if FALSE.
#' @param min_p Smallest p-value displayed. Values below this threshold are set
#'   to it for plotting. Defaults to 1e-16.
#'
#' @returns A ggplot QQ plot with a uniform-null confidence band and identity line.
#' @export
#'
#' @examples
#' make_guideseq_qq_plot_simple(stats::runif(100))
make_guideseq_qq_plot_simple <- function(p_values, rev_log_trans = TRUE, min_p = 1e-16) {
  trans <- if (rev_log_trans) revlog_trans(base = 10) else scales::reverse_trans()
  ggplot2::ggplot(data.frame(p_value = pmax(p_values, min_p)),
                  ggplot2::aes(y = p_value)) +
    stat_qq_band() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    ggplot2::scale_x_continuous(trans = trans) +
    ggplot2::scale_y_continuous(trans = trans, limits = c(1, min_p)) +
    stat_qq_points(size = 0.8, ymin = min_p, col = "black")
}
