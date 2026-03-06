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


#' Make scatterplot
#'
#' @param count_df
#' @param x_range
#' @param facet_on_chr
#' @param log_trans
#' @param point_size
#' @param col
#' @param title
#' @param target_seq
#' @param grna_spacer
#' @param pam_strand
#' @param protospacer_range
#' @param pam_range
#' @param cut_base_range
#'
#' @export
make_scatterplot <- function(count_df, x_range = NULL, facet_on_chr = FALSE, log_trans = FALSE, point_size = 1,
                             col = c("dodgerblue3", "firebrick")[1], title = NULL, target_seq = NULL,
                             grna_spacer = NULL, pam_strand = NULL, protospacer_range = NULL, pam_range = NULL,
                             cut_base_range = NULL) {
  # base plot
  p <- ggplot2::ggplot(data = count_df, mapping = ggplot2::aes(x = coord, y = count)) +
    ggplot2::geom_segment(ggplot2::aes(x = coord, xend = coord, y = if (log_trans) 1 else 0, yend = count)) +
    ggplot2::geom_point(size = point_size, col = col) +
    ggplot2::theme_bw(base_size = 10) + ggplot2::xlab("Coordinate") +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank())

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
    p <- p + ggplot2::ylab("UMI count (linear)") # + ggplot2::scale_y_continuous(expand = c(0.01, 1))
  }
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  # primary sequence info
  if (!is.null(target_seq)) {
    if (is(target_seq, "DNAString")) {
      target_seq <- strsplit(as.character(target_seq), split = "")[[1]]
    }
    if (is(grna_spacer, "DNAString")) {
      if (pam_strand == "-") grna_spacer <- Biostrings::reverseComplement(grna_spacer)
      grna_spacer <- strsplit(as.character(grna_spacer), split = "")[[1]]
    }
    max_count <- max(count_df$count)
    if (!is.null(grna_spacer)) {
      spacer_seq_y <- if (log_trans) 0.8 else -0.05 * max_count
      target_seq_y <- if (log_trans) 0.65 else -0.1 * max_count
      bottom_cut_line_y <- if (log_trans) 0.58 else -0.13 * max_count
    } else {
      target_seq_y <- if (log_trans) 0.8 else -0.05 * max_count
    }
    label_df <- data.frame(x = seq(x_range[1], x_range[2])) |>
      dplyr::mutate(base = target_seq, y = target_seq_y) |>
      dplyr::mutate(base_type = "target")

    if (!is.null(grna_spacer)) {
      label_df <- label_df |>
        dplyr::mutate(pam_site = (x >= pam_range[1] & x <= pam_range[2])) |>
        dplyr::mutate(base_type = ifelse(pam_site, "pam", "target"),
                      pam_site = NULL)
      spacer_label_df <- data.frame(x = seq(protospacer_range[1], protospacer_range[2]),
                                    base = grna_spacer, y = spacer_seq_y, base_type = "spacer")
      label_df <- rbind(label_df, spacer_label_df)
    }

    p <- p +
      ggplot2::geom_text(ggplot2::aes(x = x, y = y, label = base, col = base_type), data = label_df) +
      ggplot2::scale_color_manual(values = c("pam" = "red", "target" = "black", "spacer" = "darkorchid")) +
      ggplot2::theme(legend.position = "none")

    if (!is.null(cut_base_range)) {
      cut_spot_x <- mean(cut_base_range)
      p <- p + ggplot2::geom_segment(mapping = ggplot2::aes(x = cut_spot_x, xend = cut_spot_x,
                                                       y = if (log_trans) 1 else 0, yend = bottom_cut_line_y),
                                data = data.frame(), linetype = "dashed", col = "orange")
    }
  }

  # custom x-range?
  if (!is.null(x_range)) {
    p <- p + ggplot2::scale_x_continuous(breaks = seq(x_range[1], x_range[2], by = 1),
                                         limits = c(x_range[1], x_range[2]))
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

  # check if discoveries are present; if so, draw line
  if (any(result_df$significant_hit)) {
    rejection_thresh <- (result_df |> plyranges::filter(significant_hit))$umi_count |> min()
  } else {
    rejection_thresh <- NULL
  }

  # shift
  x_range <- x_range + 1
  y <- y + 1
  shifted_density_df$count <- shifted_density_df$count + 1

  # create untrans histogram
  p_model_untrans <- ggplot2::ggplot(data = data.frame(y = y), mapping = ggplot2::aes(x = y)) +
    ggplot2::geom_vline(xintercept = rejection_thresh, col = "blue", linetype = "dashed") +
    ggplot2::geom_histogram(binwidth = 1, col = "black", fill = "white") +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
      data = shifted_density_df,
      ggplot2::aes(x = count, y = expected),
      linewidth = 0.9,
      color = "firebrick") +
    ggplot2::ylab("Frequency") +
    ggplot2::ggtitle("Linear y-axis") + ggplot2::xlab("UMI count")

  p_model_trans <- p_model_untrans +
    ggplot2::scale_y_continuous(transform = scales::pseudo_log_trans(), breaks = c(0, 10^seq(0, 5))) +
    ggplot2::ggtitle("Log y-axis")
  out <- list(plot_untrans = p_model_untrans, plot_trans = p_model_trans)
  return(out)
}

make_manhattan_plot <- function(res_df) {
  p <- ggplot2::ggplot(data = res_df |>
                         dplyr::mutate(p_value = ifelse(p_value < 1e-300, 1e-300, p_value)),
                       mapping = ggplot2::aes(x = coord, y = p_value)) +
    ggplot2::geom_point(cex = 0.5) +
    ggplot2::scale_y_continuous(transform = revlog_trans(),
                                breaks = 10^(-seq(0, 1000, by = 100))) +
    ggplot2::facet_wrap("chr") + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::xlab("Coordinate") + ggplot2::ylab("P-value")
  return(p)
}

make_p_value_histogram <- function(clustered_res) {
  p_thresh <- if (sum(clustered_res$significant_hit) >= 1) {
    clustered_res |> dplyr::filter(significant_hit) |> dplyr::pull(p_combined) |> max()
  } else {
    NULL
  }
  p_vals <- clustered_res$p_combined
  p_model_untrans <- ggplot2::ggplot(data = data.frame(x = p_vals), mapping = ggplot2::aes(x = x)) +
    ggplot2::geom_histogram(binwidth = 0.01, col = "black", fill = "black") +
    ggplot2::theme_bw() + ggplot2::ylab("Frequency") + ggplot2::xlab("P-value") +
    ggplot2::geom_vline(xintercept = p_thresh, col = "firebrick3")
  return(p_model_untrans)
}

make_discovery_site_scatterplots_dm <- function(res_df, clustered_res_df, plot_window_size = 24, col = c("dodgerblue3", "firebrick")[1]) {
  signif_group_ids <- clustered_res_df |> dplyr::filter(significant_hit) |> dplyr::pull(group_id)
  if (length(signif_group_ids) >= 1) {
    # loop through significant group ids
    plot_list <- lapply(signif_group_ids, function(curr_group_id) {
      curr_res_sub <- res_df |> dplyr::filter(group_id == curr_group_id)
      curr_lead_base <- curr_res_sub$coord[which.max(curr_res_sub$count)]
      curr_chr <- curr_res_sub$chr[which.max(curr_res_sub$chr)]
      count_df_sub <- res_df |>
        dplyr::mutate(coord = coord - curr_lead_base) |>
        dplyr::filter(coord >= -plot_window_size/2 & coord <= plot_window_size/2 & chr == curr_chr) |>
        dplyr::select(chr, coord, count)
      plot_title <- paste0("Chr", curr_chr, ":", curr_lead_base)
      p_log <- make_scatterplot(count_df = count_df_sub,
                                x_range = c(-plot_window_size/2, plot_window_size/2),
                                facet_on_chr = FALSE, log_trans = TRUE, col = col,
                                title = plot_title)
      p_linear <- make_scatterplot(count_df = count_df_sub,
                                   x_range = c(-plot_window_size/2, plot_window_size/2),
                                   facet_on_chr = FALSE, log_trans = FALSE, col = col,
                                   title = plot_title)
      list(p_log = p_log, p_linear = p_linear)
    })
    log_plots <- lapply(X = plot_list, FUN = function(l) l[["p_log"]]) |> setNames(paste0("group:", signif_group_ids))
    linear_plots <- lapply(X = plot_list, FUN = function(l) l[["p_linear"]]) |> setNames(paste0("group:", signif_group_ids))
    out <- list(log_plots = log_plots, linear_plots = linear_plots)
  } else {
    out <- list()
  }
  return(out)
}


make_scatterplot_v2 <- function(count_df) {
  # set up labels
  chr_levels <- c(seq(1L, 22L), c("X", "Y", "M"))
  to_plot <- count_df |>
    dplyr::mutate(base_type = ifelse(on_target, "On target", "Off target"),
                  chr = gsub(x = chr, pattern = "chr", replacement = "") |>
                    factor(levels = chr_levels, labels = chr_levels))

  p <- ggplot2::ggplot(to_plot,
                       mapping = ggplot2::aes(x = coord, y = count, col = base_type)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::theme_bw() + ggplot2::xlab("Chromosome position") +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank()) +
    ggplot2::facet_grid(. ~ chr, scales = "free_x", drop = FALSE) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.spacing = grid::unit(0, "lines"),
                   strip.placement = "outside",
                   strip.background = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_text(),
                   panel.border = ggplot2::element_rect(color = "grey80", fill = NA),
                   legend.title = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   legend.position = "bottom") +
    ggplot2::ylab("UMI count") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.5)) +
    ggplot2::scale_color_manual(values = c("black", "forestgreen"))
}
