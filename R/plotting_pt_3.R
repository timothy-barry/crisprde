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
make_guideseq_qq_plot <- function(res_df, color_ground_truth = FALSE) {
  rejected_p_vals <- res_df |>
    dplyr::filter(nominated_off_target) |>
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
    ggplot2::scale_x_continuous(trans = revlog_trans(base = 10)) +
    ggplot2::scale_y_continuous(trans = revlog_trans(base = 10)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    ggplot2::geom_hline(yintercept = rejection_threshold, col = "blue", linetype = "dashed")
  if (color_ground_truth) {
    p <- p + stat_qq_points(ymin = 1e-8) +
      ggplot2::scale_color_manual(values = c("firebrick1", "black")) +
      ggplot2::labs(color = "Ground truth") +
      ggplot2::scale_size_manual(values = c("Truly edited" = 1.5, "Truly unedited" = 0.5)) +
      ggplot2::guides(size = "none")
  } else {
    p <- p + stat_qq_points(ymin = 1e-8, size = 0.5)
  }
  return(p)
}
