#' Make a weighted p-value QQ plot
#'
#' Compares observed weighted p-values with their theoretical null quantiles,
#' conditional on the supplied weights.
#'
#' For fixed weights `weights`, let `P_i` be independent Uniform(0, 1) null
#' p-values and let `Q_i = min(P_i / weights_i, 1)`. For `t < 1`, the average
#' conditional null CDF of the weighted p-values is
#' `F(t) = mean(pmin(weights * t, 1))`, with `F(1) = 1`. The plot compares the
#' ordered observed weighted p-values with `F^-1(r / m)` for ranks 1 through
#' `m`.
#'
#' Observed weighted p-values are clipped at 1 before ordering, and points equal
#' to 1 are omitted from the plot. The grey ribbon is a pointwise 95 percent
#' Monte Carlo confidence band for the smallest `n_conf_windows` order
#' statistics. It conditions on the supplied weights and assumes independent
#' Uniform(0, 1) input p-values.
#'
#' @param weighted_p_vals Weighted p-values to plot. Values above 1 are clipped
#'   at 1.
#' @param weights Fixed positive weights defining the weighted null
#'   distribution. These are typically normalized to have mean 1.
#' @param point_col Point color.
#' @param B Number of Monte Carlo samples used to construct the confidence band.
#' @param n_conf_windows Number of smallest-p-value ranks included in the
#'   confidence band.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' set.seed(4)
#' w <- sample(x = c(1, 0.5, 0.1, 0.01), replace = TRUE, size = 5000)
#' w_tilde <- w / mean(w)
#' weighted_p_vals <- runif(n = 5000, min = 0, max = 1) / w_tilde
#' qq_plot <- make_weighted_qq_plot(weighted_p_vals, w_tilde, B = 100L)
make_weighted_qq_plot <- function(weighted_p_vals, weights, point_col = "black", B = 10000L, n_conf_windows = 500L) {
  # compute the theoretical null quantiles
  null_qq_quantiles <- compute_weighted_null_qq_quantiles(weights)

  # construct data frame
  qq_df <- data.frame(rank = seq_along(null_qq_quantiles),
                      theoretical_quantile = null_qq_quantiles, # quantiles of the average null distribution
                      observed_quantile = sort(pmin(weighted_p_vals, 1))) # order statistics of the p-values

  # get confidence band df; update qq_df
  confidence_band_df <- compute_monte_carlo_confidence_band(weights = weights, B = B, n_conf_windows = n_conf_windows)
  qq_df <- dplyr::left_join(qq_df, confidence_band_df, by = "rank")

  # filter out weighted p-value 1 points
  qq_df <- qq_df |> dplyr::filter(observed_quantile < 1)

  # make plot
  p <- ggplot2::ggplot(data = qq_df,
                       mapping = ggplot2::aes(x = theoretical_quantile, y = observed_quantile)) +
    ggplot2::geom_ribbon(data = qq_df |> dplyr::filter(!is.na(lower_ci), !is.na(upper_ci)),
                         mapping = ggplot2::aes(x = theoretical_quantile,
                                                ymin = lower_ci, ymax = upper_ci),
                         inherit.aes = FALSE, fill = "grey85") +
    ggplot2::geom_abline(color = "black") +
    ggplot2::geom_point(size = 0.8, col = point_col) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(trans = revlog_trans(base = 10)) +
    ggplot2::scale_y_continuous(trans = revlog_trans(base = 10)) +
    ggplot2::labs(x = "Expected weighted-null p-value", y = "Observed weighted p-value")

  return(p)
}

compute_monte_carlo_confidence_band <- function(weights, B = 5000L, n_conf_windows = 500L, alpha = 0.05) {
  set.seed(4)
  n_windows <- length(weights)
  n_conf_windows <- min(n_conf_windows, n_windows)
  null_dist_mat <- replicate(n = B, expr = {
    null_weighted_p_vals <- runif(n = n_windows, min = 0, max = 1)
    null_weighted_weighted_p_vals <- pmin(null_weighted_p_vals/weights, 1)
    null_weighted_weighted_p_vals <- sort(null_weighted_weighted_p_vals, partial = n_conf_windows)[seq_len(n_conf_windows)]
    sort(null_weighted_weighted_p_vals)
  }, simplify = TRUE) |> t()
  emp_ci_tab <- apply(X = null_dist_mat, MARGIN = 2, FUN = function(col) {
    quantile(x = col, probs = c(alpha/2, 1 - alpha/2))
  }) |> t() |> as.data.frame()
  colnames(emp_ci_tab) <- c("lower_ci", "upper_ci")
  out <- emp_ci_tab |> dplyr::mutate(rank = seq_len(n_conf_windows))
  return(out)
}

compute_weighted_null_qq_quantiles <- function(weights) {
  m <- length(weights)
  weight_groups <- rle(sort(weights[weights > 1], decreasing = TRUE))
  saturated_n <- c(0L, cumsum(weight_groups$lengths))
  unsaturated_weight <- c(sum(weights), sum(weights) - cumsum(weight_groups$values * weight_groups$lengths))
  breakpoints <- c(0, 1 / weight_groups$values)
  cdf_breakpoints <- (saturated_n + breakpoints * unsaturated_weight) / m
  cdf_limit <- (tail(saturated_n, 1L) + tail(unsaturated_weight, 1L)) / m

  ranks <- seq_len(m)
  probabilities <- ranks / m
  inverse_cdf <- rep(1, m)
  below_jump <- probabilities < cdf_limit
  segments <- findInterval(probabilities[below_jump], cdf_breakpoints)
  inverse_cdf[below_jump] <- (ranks[below_jump] - saturated_n[segments]) /
    unsaturated_weight[segments]
  return(inverse_cdf)
}
