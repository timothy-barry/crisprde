make_scatterplot <- function(count_df, x_range = NULL, facet_on_chr = FALSE, log_trans = FALSE,
                             col = c("dodgerblue3", "firebrick")[1], title = NULL) {
  # base plot
  p <- ggplot2::ggplot(data = count_df, mapping = ggplot2::aes(x = coord, y = count)) +
    ggplot2::geom_segment(ggplot2::aes(x = coord, xend = coord, y = 1, yend = count)) +
    ggplot2::geom_point(size = 0.7, col = col) +
    ggplot2::theme_bw(base_size = 10) + ggplot2::xlab("Coordinate")

  # facet?
  if (facet_on_chr) {
    p <- p + ggplot2::facet_wrap("chr") + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                         axis.ticks.x = ggplot2::element_blank())
  }
  # log transform?
  if (log_trans) {
    p <- p + ggplot2::scale_y_continuous(trans = scales::log10_trans()) +
      ggplot2::ylab("UMI count (log)")
  } else {
    p <- p + ggplot2::ylab("UMI count (linear)")
  }
  # custom x-range?
  if (!is.null(x_range)) {
    suppressMessages(p <- p + ggplot2::xlim(x_range) + ggplot2::scale_x_continuous(breaks = seq(x_range[1], x_range[2], by = 1)))
  }
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  return(p)
}

make_discovery_site_scatterplots <- function(count_df, result_df, plot_window_size = 30, col = c("dodgerblue3", "firebrick")[1]) {
  # find the significant discoveries
  result_df_sig <- result_df |> dplyr::filter(significant_hit)
  if (nrow(result_df_sig) >= 1) {
    plot_list <- lapply(seq(1, nrow(result_df_sig)), function(i) {
      curr_lead_base <- result_df_sig$lead_base[i]
      curr_chr <- as.character(result_df_sig$chr[i])
      count_df_sub <- count_df |>
        dplyr::mutate(coord = coord - curr_lead_base) |>
        dplyr::filter(coord >= -plot_window_size/2 & coord <= plot_window_size/2 & chr == curr_chr)
      p_log <- make_scatterplot(count_df = count_df_sub,
                                x_range = c(-plot_window_size/2, plot_window_size/2),
                                facet_on_chr = FALSE, log_trans = TRUE, col = col,
                                title = paste0("Chr", curr_chr, ":", curr_lead_base))
      p_linear <- make_scatterplot(count_df = count_df_sub,
                                   x_range = c(-plot_window_size/2, plot_window_size/2),
                                   facet_on_chr = FALSE, log_trans = FALSE, col = col,
                                   title = paste0("Chr", curr_chr, ":", curr_lead_base))
      list(p_log = p_log, p_linear = p_linear)
    })
    lead_base_names <- paste0("Chr", result_df_sig$chr, ":", result_df_sig$lead_base)
    log_plots <- lapply(X = plot_list, FUN = function(l) l[["p_log"]]) |> setNames(lead_base_names)
    linear_plots <- lapply(X = plot_list, FUN = function(l) l[["p_linear"]]) |> setNames(lead_base_names)
    out <- list(log_plots = log_plots, linear_plots = linear_plots)
  } else {
    out <- list()
  }

  return(out)
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
    rejection_thresh <- result_df |> dplyr::filter(significant_hit) |> dplyr::pull(umi_count) |> min()
    p_model_untrans <- p_model_untrans + ggplot2::geom_vline(xintercept = rejection_thresh, col = "blue", linetype = "dashed")
  }
  p_model_trans <- p_model_untrans +
    ggplot2::scale_y_continuous(transform = scales::pseudo_log_trans(), breaks = c(0, 10^seq(0, 5))) +
    ggplot2::ggtitle("Log y-axis")
  out <- list(plot_untrans = p_model_untrans, plot_trans = p_model_trans)
  return(out)
}
