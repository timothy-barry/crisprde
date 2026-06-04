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


make_dunn_smyth_qq_plot <- function(y, fit, seed = 1) {
  set.seed(seed)
  y <- y - 1L
  mu <- unname(fit[["mu"]])
  theta <- unname(fit[["theta"]])

  p_lo <- ifelse(y <= 0, 0, pnbinom(q = y - 1L, mu = mu, size = theta))
  p_hi <- pnbinom(q = y, mu = mu, size = theta)
  u <- stats::runif(length(y), min = p_lo, max = p_hi)
  resid <- stats::qnorm(u)
  ord <- order(resid)

  qq_df <- tibble::tibble(
    theoretical = stats::qnorm(stats::ppoints(length(resid))),
    sample = resid[ord],
  )

  p <- ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick") +
    ggplot2::xlab("Theoretical normal quantiles") +
    ggplot2::ylab("Dunn-Smyth residual quantiles") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.06)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.06)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 12, 10, 10))
  p
}
