revlog_trans <- function(base = exp(1)) {
  trans <- function(x) {
    -log(x, base)
  }
  inv <- function(x) {
    base^(-x)
  }
  scales::trans_new(
    name = paste("revlog-", base, sep = ""),
    transform = trans,
    inverse = inv,
    breaks = scales::log_breaks(base = base),
    domain = c(1e-100, Inf)
  )
}

make_local_scatterplot <- function(count_df, x_range = NULL, log_trans = FALSE, point_size = 1,
                                   col = c("dodgerblue3", "firebrick")[1], title = NULL, target_seq = NULL,
                                   grna_spacer = NULL, pam_strand = NULL, protospacer_range = NULL, pam_range = NULL,
                                   cut_base_range = NULL) {
  count_df_plus <- count_df |> dplyr::filter(strand == "+")
  count_df_minus <- count_df |> dplyr::filter(strand == "-")

  # construct label_df
  if (is(target_seq, "DNAString")) {
    target_seq <- strsplit(as.character(target_seq), split = "")[[1]]
  }
  if (is(grna_spacer, "DNAString")) {
    if (pam_strand == "-") grna_spacer <- Biostrings::reverseComplement(grna_spacer)
    grna_spacer <- strsplit(as.character(grna_spacer), split = "")[[1]]
  }
  label_df <- data.frame(x = seq(x_range[1], x_range[2])) |>
    dplyr::mutate(base = target_seq, y = -0.1) |>
    dplyr::mutate(base_type = "target")
  if (!is.null(grna_spacer)) {
    label_df <- label_df |>
      dplyr::mutate(pam_site = (x >= pam_range[1] & x <= pam_range[2])) |>
      dplyr::mutate(base_type = ifelse(pam_site, "pam", "target"),
                    pam_site = NULL)
    spacer_label_df <- data.frame(x = seq(protospacer_range[1], protospacer_range[2]),
                                  base = grna_spacer, y = 0.1, base_type = "spacer")
    label_df <- rbind(label_df, spacer_label_df)
  }

  # base plot
  make_base_plot <- function(curr_count_df) {
    ggplot2::ggplot(data = curr_count_df, mapping = ggplot2::aes(x = coord, y = count)) +
      ggplot2::geom_segment(ggplot2::aes(x = coord, xend = coord, y = if (log_trans) 1 else 0, yend = count)) +
      ggplot2::geom_point(size = point_size, col = col) +
      ggplot2::theme_bw(base_size = 10) + ggplot2::xlab("Coordinate") +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.x  = element_blank(),
                     axis.ticks.x = element_blank(),
                     panel.border = element_blank(),
                     plot.margin = margin(0.0, 5.5, 0.0, 5.5)) +
      ggplot2::scale_x_continuous(limits = range(label_df$x))
  }
  p_plus <- make_base_plot(count_df_plus)
  p_minus <- make_base_plot(count_df_minus)

  # y-axis scale
  y_max <- max(max(count_df_plus$count), max(count_df_minus$count))
  y_limits <- if (log_trans) c(1, y_max) else c(0, y_max)
  if (log_trans) {
    p_plus <- p_plus +
      ggplot2::scale_y_continuous(trans = scales::log10_trans(),
                                  expand = c(0.01, 0),
                                  limits = y_limits) +
      ggplot2::ylab("")
    p_minus <- p_minus +
      ggplot2::scale_y_continuous(trans = revlog_trans(base = 10),
                                  expand = c(0.01, 0),
                                  limits = y_limits[2:1]) +
      ggplot2::ylab("")
  } else {
    p_plus <- p_plus + ggplot2::ylab("") +
      ggplot2::scale_y_continuous(expand = expansion(mult = 0.01, add = 0),
                                  limits = y_limits)
    p_minus <- p_minus + ggplot2::ylab("") +
      ggplot2::scale_y_continuous(trans = scales::reverse_trans(),
                                  expand = expansion(mult = 0.01, add = 0),
                                  limits = y_limits)
  }

  p_middle <- ggplot2::ggplot() +
    ggplot2::geom_text(ggplot2::aes(x = x, y = y, label = base, col = base_type), data = label_df) +
    ggplot2::scale_color_manual(values = c("pam" = "red", "target" = "black", "spacer" = "darkorchid")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x  = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.y  = element_blank(),
                   axis.ticks.y = element_blank(),
                   panel.border = element_blank(),
                   plot.margin = margin(0.0, 5.5, 0.0, 5.5),
                   legend.position = "none") +
    ylab(paste0("UMI count ", if (log_trans) "(log)" else "(linear)")) +
    scale_y_continuous(limits = c(-0.17, 0.15))

  if (!is.null(cut_base_range)) {
    cut_spot_x <- mean(cut_base_range)
    p_middle <- p_middle + ggplot2::geom_segment(mapping = ggplot2::aes(x = cut_spot_x, xend = cut_spot_x,
                                                                        y = -0.5, yend = 0.5),
                                                 data = data.frame(), col = "orange")

  }

  library(patchwork)
  p_all <- (p_plus /  p_middle / p_minus) +
    plot_layout(heights = c(1, 0.1, 1), axes = "collect") +
    plot_annotation(title = title)
}
