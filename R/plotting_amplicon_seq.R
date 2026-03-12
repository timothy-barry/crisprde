make_pilot_dispersion_plot <- function(res, label_outliers = FALSE) {
  result_df <- res$result_df
  diagnostics <- res$dispersion_diagnostics

  if (is.null(diagnostics)) {
    stop("Expected `res$dispersion_diagnostics` to be present.")
  }

  pilot_df <- result_df |>
    dplyr::filter(!is.na(pilot_rho_hat)) |>
    dplyr::mutate(amplicon_id = as.character(amplicon_id))
  final_df <- result_df |>
    dplyr::filter(!is.na(rho_hat)) |>
    dplyr::mutate(amplicon_id = as.character(amplicon_id))
  outlier_df <- result_df |>
    dplyr::filter(dispersion_outlier) |>
    dplyr::mutate(amplicon_id = as.character(amplicon_id))

  n_outliers <- nrow(outlier_df)
  outlier_colors <- if (n_outliers > 0) {
    grDevices::hcl.colors(n_outliers, palette = "Dark 3")
  } else {
    character(0)
  }
  names(outlier_colors) <- outlier_df$amplicon_id

  line_df <- data.frame(
    line_label = c("Outlier threshold", "Shared rho"),
    yintercept = c(
      diagnostics$outlier_thresh,
      diagnostics$shared_rho_hat
    ),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(!is.na(yintercept)) |>
    dplyr::mutate(
      line_label = factor(
        line_label,
        levels = c("Outlier threshold", "Median pilot rho", "Shared rho")
      )
    )

  my_plot <- ggplot2::ggplot(result_df, ggplot2::aes(x = amplicon_id)) +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(
      data = line_df,
      ggplot2::aes(yintercept = yintercept, color = line_label, linetype = line_label),
      show.legend = TRUE
    ) +
    ggplot2::geom_point(
      data = pilot_df,
      ggplot2::aes(y = pilot_rho_hat, shape = "Pilot rho"),
      color = "black",
      size = 2,
      show.legend = TRUE
    ) +
    ggplot2::geom_point(
      data = final_df,
      ggplot2::aes(y = rho_hat, shape = "Final rho"),
      color = "grey55",
      size = 3,
      stroke = 0.9,
      fill = NA,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      name = NULL,
      values = c(
        "Outlier threshold" = "firebrick",
        "Median pilot rho" = "dodgerblue2",
        "Shared rho" = "goldenrod"
      )
    ) +
    ggplot2::scale_shape_manual(
      name = NULL,
      breaks = c("Pilot rho", "Final rho"),
      values = c("Pilot rho" = 16, "Final rho" = 1)
    ) +
    ggplot2::scale_linetype_manual(
      name = NULL,
      values = c(
        "Outlier threshold" = "dashed",
        "Median pilot rho" = "solid",
        "Shared rho" = "solid"
      )
    ) +
    ggplot2::scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 5e-5, base = 10),
      breaks = c(0, 1e-4, 1e-3, 1e-2, 1e-1),
      expand = ggplot2::expansion(mult = c(0.08, 0.18))
    ) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 1)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::xlab("Amplicon") +
    ggplot2::ylab(expression(hat(rho))) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        order = 1,
        override.aes = list(shape = NA, linewidth = 0.8)
      ),
      shape = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          shape = c(16, 1),
          color = c("black", "grey55"),
          fill = c("black", NA),
          size = c(2, 3),
          stroke = c(0, 0.9),
          linewidth = 0
        )
      ),
      linetype = "none",
      fill = ggplot2::guide_legend(
        order = 3,
        override.aes = list(shape = 21, size = 3, color = "black")
      )
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.margin = ggplot2::margin(10, 18, 12, 12)
    )

  if (n_outliers > 0) {
    if (label_outliers) {
      my_plot <- my_plot +
        ggplot2::geom_point(
          data = outlier_df,
          ggplot2::aes(y = pilot_rho_hat, fill = amplicon_id),
          shape = 21,
          color = "black",
          size = 2.5,
          show.legend = TRUE
        ) +
        ggplot2::scale_fill_manual(name = "Outlier amplicons", values = outlier_colors)
    } else {
      my_plot <- my_plot +
        ggplot2::geom_point(
          data = outlier_df,
          ggplot2::aes(y = pilot_rho_hat),
          shape = 21,
          fill = "black",
          color = "black",
          size = 2.5,
          show.legend = FALSE
        )
    }
  }

  my_plot
}

make_amplicon_seq_ci_plot <- function(result_df, ylim = NULL, point_size = 1.2,
                                      color_by_significance = FALSE,
                                      highlight_upper_ci_exceeds = FALSE,
                                      upper_ci_threshold = 0.001) {
  result_df <- result_df |>
    dplyr::mutate(
      Significant = significant,
      upper_ci_flag = theta_hat_upper_ci > upper_ci_threshold
    )

  if (color_by_significance && highlight_upper_ci_exceeds) {
    my_plot <- ggplot2::ggplot(
      data = result_df,
      mapping = ggplot2::aes(
        x = amplicon_id,
        y = theta_hat_clipped,
        color = Significant
      )
    ) +
      ggplot2::geom_point(size = point_size) + ggplot2::theme_bw() +
      ggplot2::theme_bw() +
      ggplot2::geom_errorbar(mapping = ggplot2::aes(ymin = theta_hat_lower_ci,
                                                    ymax = theta_hat_upper_ci, width = 0,
                                                    color = Significant)) +
      ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "dodgerblue2"))
  } else if (color_by_significance) {
    my_plot <- ggplot2::ggplot(
      data = result_df,
      mapping = ggplot2::aes(x = amplicon_id, y = theta_hat_clipped, color = Significant)
    ) +
      ggplot2::geom_point(size = point_size) + ggplot2::theme_bw() +
      ggplot2::theme_bw() +
      ggplot2::geom_errorbar(mapping = ggplot2::aes(ymin = theta_hat_lower_ci,
                                                    ymax = theta_hat_upper_ci, width = 0)) +
      ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "dodgerblue2"))
  } else if (highlight_upper_ci_exceeds) {
    my_plot <- ggplot2::ggplot(
      data = result_df,
      mapping = ggplot2::aes(x = amplicon_id, y = theta_hat_clipped)
    ) +
      ggplot2::geom_point(size = point_size) + ggplot2::theme_bw() +
      ggplot2::theme_bw() +
      ggplot2::geom_errorbar(mapping = ggplot2::aes(ymin = theta_hat_lower_ci,
                                                    ymax = theta_hat_upper_ci, width = 0))
  } else {
    my_plot <- ggplot2::ggplot(data = result_df,
                               mapping = ggplot2::aes(x = amplicon_id, y = theta_hat_clipped)) +
      ggplot2::geom_point(size = point_size) + ggplot2::theme_bw() +
      ggplot2::theme_bw() +
      ggplot2::geom_errorbar(mapping = ggplot2::aes(ymin = theta_hat_lower_ci,
                                                    ymax = theta_hat_upper_ci, width = 0))
  }

  my_plot <- my_plot +
    ggplot2::scale_y_continuous(labels = scales::label_percent(), limits = ylim) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::xlab("Amplicon") + ggplot2::ylab("Estimated editing rate")

  if (highlight_upper_ci_exceeds) {
    my_plot <- my_plot +
      ggplot2::geom_point(
        data = result_df |>
          dplyr::filter(upper_ci_flag),
        ggplot2::aes(x = amplicon_id, y = theta_hat_clipped, shape = "Upper CI exceeds safety threshold"),
        inherit.aes = FALSE,
        size = point_size + 1.1,
        color = "black"
      ) +
      ggplot2::scale_shape_manual(
        name = NULL,
        values = c("Upper CI exceeds safety threshold" = 17)
      )
  }

  my_plot
}

make_amplicon_seq_p_value_plot <- function(result_df, min_p_value = 1e-250) {
  ggplot2::ggplot(data = result_df |>
                    dplyr::mutate(Significant = significant,
                                  p_value = pmax(p_value, min_p_value)),
                  mapping = ggplot2::aes(x = amplicon_id, y = p_value, col = Significant)) +
    ggplot2::geom_point() + ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                   legend.position = "bottom") +
    ggplot2::scale_color_manual(values = c("black", "dodgerblue2")) +
    ggplot2::xlab("Amplicon") +
    ggplot2::scale_y_continuous(trans = sceptre:::revlog_trans(base = 10)) +
    ggplot2::ylab("p-value")
}

convert_data_list_into_mutation_frac_df <- function(data_list) {
  compute_mutation_frac <- function(n_mat, k_mat) {
    as.data.frame(k_mat/n_mat) |>
      tidyr::pivot_longer(cols = tidyr::everything(),
                          names_to = "amplicon", values_to = "mutation_frac") |>
      dplyr::arrange(amplicon)
  }
  mutation_frac_df <- rbind(compute_mutation_frac(data_list$n_mat_trt, data_list$k_mat_trt) |>
                              dplyr::mutate(Condition = "Treated"),
                            compute_mutation_frac(data_list$n_mat_cntrl, data_list$k_mat_cntrl) |>
                              dplyr::mutate(Condition = "Control"))

  n_vect <- colSums(data_list$n_mat_trt) + colSums(data_list$n_mat_cntrl)
  n_reads_df <- data.frame(amplicon = names(n_vect), n_reads = setNames(n_vect, NULL))

  return(list(mutation_frac_df = mutation_frac_df, n_reads_df = n_reads_df))
}

make_mutation_frac_plot <- function(to_plot) {
  ggplot2::ggplot(data = to_plot,
                  mapping = ggplot2::aes(x = Condition, y = mutation_frac, col = Condition)) +
    ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::xlab("Sample") + ggplot2::ylab("Mutation fraction") +
    ggplot2::facet_wrap(. ~ amplicon) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 8), legend.position = "bottom") +
    ggplot2::scale_color_manual(values = c("firebrick1", "dodgerblue1"))
}

make_n_reads_plot <- function(to_plot) {
  ggplot2::ggplot(data = to_plot,
                  mapping = ggplot2::aes(x = amplicon, y = n_reads)) +
    ggplot2::geom_point() + ggplot2::theme_bw()
}
