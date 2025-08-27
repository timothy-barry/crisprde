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
    rejection_thresh <- result_df |> dplyr::filter(significant_hit) |> dplyr::pull(umi_count) |> min()
    p_model_untrans <- p_model_untrans + ggplot2::geom_vline(xintercept = rejection_thresh, col = "blue", linetype = "dashed")
  }
  p_model_trans <- p_model_untrans +
    ggplot2::scale_y_continuous(transform = scales::pseudo_log_trans(), breaks = c(0, 10^seq(0, 5))) +
    ggplot2::ggtitle("Log y-axis")
  out <- list(plot_untrans = p_model_untrans, plot_trans = p_model_trans)
  return(out)
}
