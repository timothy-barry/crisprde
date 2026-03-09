make_pilot_dispersion_plot <- function(res) {
  result_df <- res$result_df
  diagnostics <- res$dispersion_diagnostics

  if (is.null(diagnostics)) {
    stop("Expected `res$dispersion_diagnostics` to be present.")
  }

  line_df <- data.frame(
    line_label = c("Outlier threshold", "Median pilot rho", "Shared rho"),
    yintercept = c(
      diagnostics$outlier_thresh,
      diagnostics$median_pilot_rho,
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

  ggplot2::ggplot(result_df, ggplot2::aes(x = amplicon_id, y = pilot_rho_hat)) +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(
      data = line_df,
      ggplot2::aes(yintercept = yintercept, color = line_label, linetype = line_label),
      show.legend = TRUE
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(
      name = NULL,
      values = c(
        "Outlier threshold" = "firebrick",
        "Median pilot rho" = "dodgerblue2",
        "Shared rho" = "goldenrod"
      )
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
      breaks = c(0, 1e-4, 1e-3, 1e-2, 1e-1)
    ) +
    ggplot2::xlab("Amplicon") +
    ggplot2::ylab(expression(hat(rho))) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
}

make_amplicon_seq_ci_plot <- function(result_df, ylim = NULL) {
  my_plot <- ggplot2::ggplot(data = result_df,
                             mapping = ggplot2::aes(x = amplicon_id, y = theta_hat_clipped)) +
    ggplot2::geom_point() + ggplot2::theme_bw() +
    ggplot2::theme_bw() +
    ggplot2::geom_errorbar(mapping = ggplot2::aes(ymin = theta_hat_lower_ci,
                                                  ymax = theta_hat_upper_ci, width = 0)) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(), limits = ylim) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::xlab("Amplicon") + ggplot2::ylab("Estimated editing rate")
}

make_amplicon_seq_p_value_plot <- function(result_df, min_p_value = 1e-250) {
  ggplot2::ggplot(data = result_df |>
                    dplyr::mutate(Significant = significant,
                                  p_value = pmax(p_value, min_p_value)),
                  mapping = ggplot2::aes(x = amplicon_id, y = p_value, col = Significant)) +
    ggplot2::geom_point() + ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::scale_color_manual(values = c("black", "dodgerblue2")) +
    ggplot2::xlab("Amplicon") + ggplot2::scale_y_continuous(trans = sceptre:::revlog_trans(base = 10)) +
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
