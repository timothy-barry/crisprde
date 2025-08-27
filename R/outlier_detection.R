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
#' @param multiplicity_alpha
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
find_guideseq_edit_sites <- function(count_df, window_size = 1000, multiplicity_adjustment = "BH", model = "nb",
                                     robust_mle = TRUE, multiplicity_alpha = 0.1, c_tukey_beta = 10, c_tukey_sigma = 10,
                                     chrs_to_keep = seq(1L, 22L)) {
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
  count_vec <- tapply(GenomicRanges::mcols(gr_reads)$count[S4Vectors::subjectHits(hits)], S4Vectors::queryHits(hits), sum)
  gr_bins$read_sum <- 0
  gr_bins$read_sum[as.integer(names(count_vec))] <- count_vec

  # 2.d keep only bins whose read count is greater than or equal to zero
  gr_bins_sub <- gr_bins[gr_bins$read_sum >= 1]
  v <- gr_bins_sub$read_sum

  # 3. shift v, then fit the model
  y <- v - 1L
  if (robust_mle && model == "nb") {
    fit <- fit_rob_nb_univariate(y = y, c.tukey.beta = c_tukey_beta, c.tukey.sigma = c_tukey_sigma)
    theta_hat <- fit[["theta"]]
  } else if (!robust_mle && model == "nb") {
    fit <- fit_nb_univariate(y = y)
    theta_hat <- fit[["theta"]]
  } else if (!robust_mle && model == "poisson") {
      fit <- list(mu = mean(y))
  } else {
    stop("Model not recognized.")
  }
  mu_hat <- fit[["mu"]]

  # 4. compute p-value for each observation
  if (model == "nb") {
    p_vals <- compute_exact_p_value_nb(mu_hat, theta_hat, y)
  } else if (model == "poisson") {
    p_vals <- compute_exact_p_value_poisson(mu_hat, y)
  } else {
    stop("Model not recognized.")
  }

  # 5. apply multiplicity adjustment
  p_adj <- p.adjust(p = p_vals, method = multiplicity_adjustment)
  rejected <- (p_adj < multiplicity_alpha)

  # 6. prepare result df
  result_df <- data.frame(chr = gr_bins_sub@seqnames,
                          range_start = gr_bins_sub@ranges@start,
                          range_end = gr_bins_sub@ranges@start + gr_bins_sub@ranges@width,
                          umi_count = gr_bins_sub$read_sum, p_value = p_vals,
                          p_adj = p_adj, significant_hit = rejected)

  # 7. create plot
  plot_list <- make_histogram_plot(y, fit, result_df, model)

  # 8. prepare output
  out <- list(result_df = result_df, fitted_model = fit,
              plot_untrans = plot_list$plot_untrans,
              plot_trans = plot_list$plot_trans)
  return(out)
}

compute_exact_p_value_poisson <- function(mu, y) {
  p_vals <- stats::ppois(q = y - 1, lambda = mu, lower.tail = FALSE)
  return(p_vals)
}

compute_exact_p_value_nb <- function(mu, theta, y) {
  p_vals <- stats::pnbinom(q = y - 1, size = theta, mu = mu, lower.tail = FALSE)
  return(p_vals)
}

compute_lrt_p_value_nb <- function(mu, theta, y) {
  xi_squared <- 2 * (dnbinom(x = y, size = theta, mu = y, log = TRUE) - dnbinom(x = y, size = theta, mu = mu, log = TRUE))
  p_vals <- ifelse(y <= mu, 1, 0.5 * pchisq(xi_squared, df = 1, lower.tail = FALSE))
  return(p_vals)
}
