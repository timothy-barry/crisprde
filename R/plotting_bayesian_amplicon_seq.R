make_bayesian_density_plot <- function(density_df, parameter, title_text, line_color,
                                       x_limits = c(0, 1), vline_x = NULL,
                                       vline_color = "red", vline_linetype = "solid") {
  x_label <- if (parameter == "theta") expression(theta) else expression(pi)

  p <- ggplot2::ggplot(density_df, ggplot2::aes(x = value, y = density)) +
    ggplot2::geom_line(color = line_color) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 100)) +
    ggplot2::coord_cartesian(xlim = x_limits) +
    ggplot2::xlab(x_label) +
    ggplot2::ylab("Density") +
    ggplot2::ggtitle(title_text)

  if (!is.null(vline_x)) {
    p <- p + ggplot2::geom_vline(xintercept = vline_x, color = vline_color, linetype = vline_linetype)
  }

  p
}


#' Plot a Beta prior density for Bayesian amplicon-seq analysis
#'
#' Draw the probability density function of the Beta prior used for either the
#' background mutation rate `pi` or the editing rate `theta`.
#'
#' @param parameter Which parameter prior to plot, either `"theta"` or `"pi"`.
#' @param alpha First shape parameter of the Beta prior.
#' @param beta Second shape parameter of the Beta prior.
#' @param n_grid Number of grid points used to evaluate the density on `(0, 1)`.
#' @param line_color Color of the density curve.
#' @param xmin Lower endpoint used to evaluate the density grid.
#' @param xmax Upper endpoint used to evaluate the density grid.
#'
#' @returns A `ggplot2` object showing the requested Beta prior density.
#' @export
#'
#' @examples
#' make_prior_density_plot(parameter = "theta", alpha = 1, beta = 1)
#' make_prior_density_plot(parameter = "pi", alpha = 2, beta = 50, line_color = "firebrick")
make_prior_density_plot <- function(alpha, beta, parameter = c("theta", "pi"),
                                    n_grid = 1000L, line_color = "dodgerblue3",
                                    xmin = .Machine$double.eps, xmax = 1 - .Machine$double.eps) {
  parameter <- match.arg(parameter)
  n_grid <- as.integer(n_grid)
  x_grid <- seq(from = xmin, to = xmax, length.out = n_grid)
  density_df <- data.frame(
    value = x_grid,
    density = stats::dbeta(x_grid, shape1 = alpha, shape2 = beta)
  )
  title_text <- paste0(
    "Prior density of ", parameter, " (alpha = ",
    signif(alpha, 4), ", beta = ", signif(beta, 4), ")"
  )

  make_bayesian_density_plot(
    density_df = density_df,
    parameter = parameter,
    title_text = title_text,
    line_color = line_color,
    x_limits = c(xmin, xmax)
  )
}


#' Plot the posterior density from Bayesian amplicon-seq analysis
#'
#' Plot the marginal posterior density of `theta` or `pi` for one amplicon.
#'
#' @param posterior_density_df Combined posterior density table from
#'   [run_bayesian_amplicon_seq_analysis()].
#' @param amplicon_id_to_plot Amplicon identifier to plot.
#' @param parameter Which posterior to plot, either `"theta"` or `"pi"`.
#' @param line_color Color of the density curve.
#' @param x_limits Length-2 numeric vector giving the x-axis plotting window.
#'
#' @returns A `ggplot2` object showing the requested posterior density.
#' @export
#'
#' @examples
#' posterior_density_df <- data.frame(
#'   theta = seq(0, 1, length.out = 200),
#'   density = stats::dbeta(seq(0, 1, length.out = 200), 15, 585),
#'   amplicon_id = "amplicon_1"
#' )
#' make_posterior_density_plot(
#'   posterior_density_df = posterior_density_df,
#'   amplicon_id_to_plot = "amplicon_1",
#'   parameter = "theta",
#'   x_limits = c(0, 0.1)
#' )
make_posterior_density_plot <- function(posterior_density_df, amplicon_id_to_plot, parameter,
                                        line_color = "darkorange3", x_limits = c(0, 1)) {
  density_df <- posterior_density_df |> dplyr::filter(amplicon_id == amplicon_id_to_plot)
  density_df$value <- density_df[[parameter]]
  title_text <- paste0("Posterior density of ", parameter, " (", amplicon_id_to_plot, ")")

  make_bayesian_density_plot(
    density_df = density_df,
    parameter = parameter,
    title_text = title_text,
    line_color = line_color,
    x_limits = x_limits
  )
}
