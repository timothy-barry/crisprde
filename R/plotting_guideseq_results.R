#' Make GUIDE-seq QQ plot
#'
#' @param res_df A data frame with raw p-values in `p_value`. Optional logical
#'   `nominated_window` marks sites for the threshold line and discovery labels;
#'   `window_label` supplies labels when `annotate_discoveries = TRUE`.
#'   Weighted p-values in `p_value_weighted` are not used; plot those with
#'   [make_weighted_qq_plot()].
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
  nominated_rows <- which(res_df[["nominated_window"]] %in% TRUE)
  rejected_p_vals <- res_df$p_value[nominated_rows]
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
  if (annotate_discoveries && length(nominated_rows) > 0L) {
    label_df <- res_df |> dplyr::mutate(
      p_value_plot = pmax(p_value, min_p),
      qq_x = stats::qunif(stats::ppoints(dplyr::n())[rank(p_value_plot, ties.method = "first")]),
      label_y = pmin(1, p_value_plot * 1.5)
    ) |>
      dplyr::slice(nominated_rows)
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


#' Plot fitted GUIDE-seq count distributions
#'
#' Plots a histogram of the counts, disaggregated by condition and channel, and
#' superimposes the fitted count curve.
#'
#' @param Y_mat_trt Treated count matrix with channels in rows and windows in columns.
#' @param Y_mat_cntrl Control count matrix with channels in rows and windows in columns.
#' @param hyperparam_res Output of [tune_hyperparameters()], including
#'   `count_model` (a model label or named channel-to-model vector) and the fitted parameters.
#' @param c Robustness parameter to display. If NULL, use each condition's selected fit.
#' @param binwidth Histogram bin width in UMIs; defaults to 10.
#'
#' @details Only positive counts are plotted. Curves show
#' expected bin counts under the fitted shifted NB or shifted NB mixture.
#' @returns A ggplot with channels in rows, control and treated conditions in
#'   columns, and a pseudo-log y-axis (`sigma = 1`, `base = 10`).
#' @export
make_fitted_count_distribution_plots <- function(Y_mat_trt, Y_mat_cntrl, hyperparam_res, c = NULL, binwidth = 10) {
  Y_mats <- list(cntrl = Y_mat_cntrl, trt = Y_mat_trt)
  count_model <- hyperparam_res$count_model
  plot_df <- lapply(names(Y_mats), function(condition) {
    channel_models <- get_channel_count_models(count_model, rownames(Y_mats[[condition]]))
    parameters <- if (is.null(c)) {
      hyperparam_res[[paste0("selected_", condition, "_run")]]$ests_list$mu_theta_hat_mat
    } else {
      hyperparam_res$nb_model_fits[[condition]][[as.character(c)]]
    }
    # Cover the last bin in full across the condition's shared x-axis.
    upper_boundary <- 0.5 + (floor((max(Y_mats[[condition]]) - 0.5) / binwidth) + 1) * binwidth
    count_grid <- seq_len(ceiling(upper_boundary) - 1L)
    lapply(rownames(Y_mats[[condition]]), function(channel) {
      y <- Y_mats[[condition]][channel, ]
      y <- y[y >= 1]
      par <- parameters[channel, ]
      pmf <- if (channel_models[[channel]] == "nb") {
        dnbinom(count_grid - 1L, mu = par[["mu"]], size = par[["theta"]])
      } else {
        (1 - par[["pi"]]) * dnbinom(count_grid - 1L,
          mu = par[["mu_lower"]], size = par[["theta_lower"]]) +
          par[["pi"]] * dnbinom(count_grid - 1L,
            mu = par[["mu_upper"]], size = par[["theta_upper"]])
      }
      dplyr::bind_rows(
        data.frame(condition, channel, umi_count = y, weight = 1, distribution = "observed"),
        data.frame(condition, channel, umi_count = count_grid,
                   weight = length(y) * pmf, distribution = "fitted")
      )
    }) |> dplyr::bind_rows()
  }) |> dplyr::bind_rows()
  plot_df$channel <- factor(plot_df$channel, levels = rownames(Y_mat_trt))
  plot_df$condition <- factor(plot_df$condition, levels = c("cntrl", "trt"),
    labels = c("Control", "Treated"))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = umi_count, weight = weight)) +
    ggplot2::geom_histogram(data = dplyr::filter(plot_df, distribution == "observed"),
      binwidth = binwidth, boundary = 0.5, closed = "left",
      fill = "grey70", colour = "black", linewidth = 0.15) +
    ggplot2::stat_bin(data = dplyr::filter(plot_df, distribution == "fitted"),
      geom = "line", position = "identity", binwidth = binwidth, boundary = 0.5, closed = "left",
      colour = "#0072B2", linewidth = 0.8) +
    ggplot2::facet_grid(channel ~ condition, scales = "free") +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10),
      breaks = 10^(0:5), labels = scales::label_comma()) +
    ggplot2::labs(x = "Window-wise UMI count", y = "Number of windows (pseudo-log scale)") +
    ggplot2::theme_bw()
  return(p)
}



#' Plot fitted GUIDE-seq occupancy distributions
#'
#' Compare observed and fitted probabilities for every nonempty occupancy pattern.
#'
#' @param Y_mat_trt Treated count matrix with named channels in rows and nonempty
#'   windows in columns.
#' @param Y_mat_cntrl Control count matrix with the same named channels and
#'   nonempty windows in columns.
#' @param hyperparam_res Output of [tune_hyperparameters()] with fitted occupancy
#'   models in `occupancy_model_fits` for both conditions.
#'
#' @details Uses the saved occupancy assignments and fitted probabilities.
#'   Bars show observed proportions, with raw window counts above them; blue
#'   points and lines show fitted probabilities. Patterns are ordered by the
#'   number of occupied channels, identified in the matrix below the plot.
#'   Zero observed proportions are marked at the display floor of 1e-9.
#' @returns A patchwork plot with control on the left, treated on the right,
#'   and a shared log10 probability scale.
#' @export
make_fitted_occupancy_distribution_plots <- function(Y_mat_trt, Y_mat_cntrl, hyperparam_res) {
  channels <- rownames(Y_mat_trt)
  n_channels <- nrow(Y_mat_trt)
  Y_mats <- list(cntrl = Y_mat_cntrl, trt = Y_mat_trt)
  n_windows <- vapply(Y_mats, ncol, integer(1))
  Omega <- hyperparam_res$occupancy_model_fits$trt$Omega
  patterns <- data.frame(pattern = apply(Omega, 1, paste0, collapse = ""),
                         occupied_channels = rowSums(Omega)) |>
    dplyr::arrange(occupied_channels, dplyr::desc(pattern)) |>
    dplyr::mutate(position = dplyr::row_number())
  dividers <- which(diff(patterns$occupied_channels) != 0) + 0.5
  x_limits <- c(0.5, nrow(patterns) + 0.5)

  plot_df <- lapply(names(Y_mats), function(condition) {
    fit <- hyperparam_res$occupancy_model_fits[[condition]]
    fitted_patterns <- apply(fit$Omega[, match(channels, names(fit$pi_hat)), drop = FALSE], MARGIN = 1,
                             FUN =  paste0, collapse = "")
    pattern_order <- match(patterns$pattern, fitted_patterns)
    observed_count <- tabulate(fit$occupancy_pattern_map, nbins = nrow(fit$Omega))[pattern_order]
    data.frame(condition, patterns, observed_count,
               observed = observed_count / n_windows[[condition]],
               expected = fit$tbp_pattern_df$pmf[pattern_order])
  }) |> dplyr::bind_rows()
  plot_df$condition <- factor(plot_df$condition, levels = names(Y_mats),
                              labels = c("Control", "Treated"))

  probability_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = position)) +
    ggplot2::geom_vline(xintercept = dividers, colour = "grey85", linetype = "dashed") +
    ggplot2::geom_rect(data = dplyr::filter(plot_df, observed > 0),
      ggplot2::aes(xmin = position - 0.35, xmax = position + 0.35,
                   ymin = 1e-9, ymax = observed, fill = "Observed"),
      colour = "grey35", linewidth = 0.2) +
    ggplot2::geom_point(data = dplyr::filter(plot_df, observed == 0),
      ggplot2::aes(y = 1e-9), shape = 4, size = 2, colour = "grey35") +
    ggplot2::geom_line(ggplot2::aes(y = expected, colour = "Fitted occupancy model"), linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(y = expected, colour = "Fitted occupancy model"), size = 2.4) +
    ggplot2::geom_text(ggplot2::aes(y = pmax(observed, 1e-9), label = scales::comma(observed_count)),
      size = 2.5, vjust = -1.6) +
    ggplot2::facet_grid(. ~ condition) +
    ggplot2::scale_fill_manual(values = c(Observed = "grey75"), name = NULL) +
    ggplot2::scale_colour_manual(values = c("Fitted occupancy model" = "#0072B2"), name = NULL) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    ggplot2::scale_y_log10(limits = c(1e-9, 1), breaks = 10^seq(-9, 0),
      labels = parse(text = paste0("10^", seq(-9, 0)))) +
    ggplot2::labs(x = NULL, y = "Probability (log10 scale)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95"))

  matrix_df <- expand.grid(position = patterns$position, channel = seq_len(n_channels),
                           condition = factor(levels(plot_df$condition), levels = levels(plot_df$condition))) |>
    dplyr::mutate(occupied = substr(patterns$pattern[position], channel, channel) == "1")
  connections <- matrix_df |>
    dplyr::filter(occupied) |>
    dplyr::group_by(condition, position) |>
    dplyr::summarise(first = min(channel), last = max(channel), .groups = "drop")
  pattern_matrix <- ggplot2::ggplot(matrix_df, ggplot2::aes(x = position, y = channel)) +
    ggplot2::geom_vline(xintercept = dividers, colour = "grey85", linetype = "dashed") +
    ggplot2::geom_point(colour = "grey88", size = 2.8) +
    ggplot2::geom_segment(data = connections,
      ggplot2::aes(xend = position, y = first, yend = last), colour = "grey25", linewidth = 0.6) +
    ggplot2::geom_point(data = dplyr::filter(matrix_df, occupied), colour = "grey25", size = 2.8) +
    ggplot2::facet_grid(. ~ condition) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    ggplot2::scale_y_reverse(breaks = seq_len(n_channels), labels = channels,
      limits = c(n_channels + 0.5, 0.5), expand = c(0, 0)) +
    ggplot2::labs(x = "Occupancy pattern", y = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(), strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank())

  p <- patchwork::wrap_plots(probability_plot, pattern_matrix, ncol = 1,
                            heights = c(4, 0.2875 * n_channels), guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  return(p)
}
