#' Find GUIDE-seq edit sites
#'
#' TO DO:
#' 1. optimize c_tukey_beta and c_tukey_sigma for our setting
#' 2. consider binning-free approaches
#' 3. Multiple samples and covariates
#'
#' @param c_tukey_beta
#' @param c_tukey_sigma
#' @param count_df
#' @param window_size
#' @param p_val_calculation_method
#' @param multiplicity_adjustment
#' @param multiplicity_adjustment_threshold
#' @param chrs_to_keep
#'
#' @returns
#' @export
#'
#' @examples
#' umi_tab_fp <- "/Users/timbarry/research_offsite/external/crispr-quant/guideseq/count_tables/293T-SpRY-Cas9-dsODN-only-GSPneg-S81_S89_L001_count_table.rds"
#' count_df <- readRDS(umi_tab_fp)
#' res_cntrl <- find_guideseq_edit_sites(count_df)
#'
#' umi_tab_fp <- "/Users/timbarry/research_offsite/external/crispr-quant/guideseq/count_tables/293T-SpRY-Cas9-1620-GSPneg-S79_S87_L001_count_table.rds"
#' count_df <- readRDS(umi_tab_fp)
#' res_trt <- find_guideseq_edit_sites(count_df)
#'
find_guideseq_edit_sites <- function(count_df, window_size = 1000, p_val_calculation_method = "exact", multiplicity_adjustment = "BH", robust_mle = TRUE,
                                     multiplicity_adjustment_threshold = 0.1, c_tukey_beta = 10, c_tukey_sigma = 10, chrs_to_keep = seq(1L, 22L)) {
  # 1. subset count_df, keeping only chromosomes in chrs_to_keep
  count_df <- count_df |> dplyr::filter(chr %in% chrs_to_keep)

  # 2. partition genome into nonoverlapping windows; compute the umi count in each window
  # 2.a make granges object
  gr_reads <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(count_df$chr),
                                     ranges = IRanges::IRanges(start = count_df$coord, width = 1), strand = "*",
                                     count = count_df$count)

  # 2.b create genome tile
  window_size <- window_size
  chrom_ranges <- count_df |> dplyr::group_by(chr) |>
    dplyr::summarize(min_coord = min(coord), max_coord = max(coord))
  gr_bins_list <- lapply(X = seq(1, nrow(chrom_ranges)), function(i) {
    curr_min_coord <- chrom_ranges[[i, "min_coord"]]
    curr_max_coord <- chrom_ranges[[i, "max_coord"]]
    chr_name <- as.character(chrom_ranges[[i, "chr"]])
    bin_points <- seq(from = curr_min_coord, to = curr_max_coord, by = window_size)
    gr_bins <- GenomicRanges::GRanges(
      seqnames = chr_name,
      ranges = IRanges::IRanges(start = bin_points, width = window_size),
      seqinfo = Seqinfo::Seqinfo(seqnames = as.character(chrom_ranges$chr))
    )
  })
  gr_bins <- do.call(c, gr_bins_list)

  # 2.c count overlaps between gr_bins and gr to get count distribution
  hits <- GenomicRanges::findOverlaps(gr_bins, gr_reads)
  count_vec <- tapply(mcols(gr_reads)$count[subjectHits(hits)], S4Vectors::queryHits(hits), sum)
  gr_bins$read_sum <- 0
  gr_bins$read_sum[as.integer(names(count_vec))] <- count_vec

  # 2.d keep only bins whose read count is greater than or equal to zero
  gr_bins_sub <- gr_bins[gr_bins$read_sum >= 1]
  v <- gr_bins_sub$read_sum

  # 3. shift v, then fit the robust m-estimator
  y <- v - 1L
  if (robust_mle) {
    fit <- fit_rob_nb_univariate(y = y, c.tukey.beta = c_tukey_beta, c.tukey.sigma = c_tukey_sigma)
  } else {
    fit <- fit_nb_univariate(y = y)
  }
  mu_hat <- fit[["mu"]]
  theta_hat <- fit[["theta"]]

  # 4. compute p-value for each observation
  if (p_val_calculation_method == "exact") {
    p_vals <- compute_exact_p_value(mu_hat, theta_hat, y)
  } else if (p_val_calculation_method == "lrt") {
    p_vals <- compute_lrt_p_value(mu_hat, theta_hat, y)
  } else {
    stop("`p_val_calculation_method` not recognized.")
  }

  # 5. apply multiplicity adjustment
  p_adj <- p.adjust(p = p_vals, method = multiplicity_adjustment)
  rejected <- p_adj < multiplicity_adjustment_threshold

  # 6. prepare result df
  result_df <- data.frame(chr = gr_bins_sub@seqnames,
                          range_start = gr_bins_sub@ranges@start,
                          range_end = gr_bins_sub@ranges@start + gr_bins_sub@ranges@width,
                          umi_count = gr_bins_sub$read_sum, p_value = p_vals,
                          p_adj = p_adj, significant_hit = rejected)

  # 7. create plot


  # 8. prepare output
  out <- list(result_df = result_df, fit = fit)
}

compute_exact_p_value <- function(mu, theta, y) {
  p_vals <- stats::pnbinom(q = y - 1, size = theta, mu = mu, lower.tail = FALSE)
  return(p_vals)
}

compute_lrt_p_value <- function(mu, theta, y) {
  xi_squared <- 2 * (dnbinom(x = y, size = theta, mu = y, log = TRUE) - dnbinom(x = y, size = theta, mu = mu, log = TRUE))
  p_vals <- ifelse(y <= mu, 1, 0.5 * pchisq(xi_squared, df = 1, lower.tail = FALSE))
  return(p_vals)
}

make_plot <- function(y, fit, result_df) {
  # compute fitted density
  x_range <- seq(0, max(y))
  shifted_nb_df <- data.frame(
    count = x_range,
    expected = dnbinom(x = x_range, size = fit[["theta"]], mu = fit[["mu"]]) * length(y)
  )

  # shift
  x_range <- x_range + 1
  y <- y + 1
  shifted_nb_df$count <- shifted_nb_df$count + 1

  # plot
  p_model_untrans <- ggplot2::ggplot(data = data.frame(y = y), mapping = ggplot2::aes(x = y)) +
    ggplot2::geom_histogram(binwidth = 1, col = "black", fill = "white") +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
      data = shifted_nb_df,
      ggplot2::aes(x = count, y = expected),
      linewidth = 0.9,
      color = "firebrick"
    ) + ggplot2::ylab("Frequency")

  p_model_trans <- p_model_untrans + ggplot2::scale_y_continuous(transform = scales::pseudo_log_trans(), breaks = c(0, 10^seq(0, 5)))
  p_model_untrans <- p_model_untrans
  p_model_trans <- p_model_trans

  p_all <- cowplot::plot_grid(p_model_untrans, p_model_trans, nrow = 1)
  p_all
}
