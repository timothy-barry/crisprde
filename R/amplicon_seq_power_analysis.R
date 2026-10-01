#' Design n and r in an amplicon-seq experiment
#'
#' Given key parameters of an amplicon-seq study, design n (number of reads) and r (number of independent samples)
#'
#' @param pi_cntrl mutation rate in control condition
#' @param theta editing rate
#' @param rho overdispersion parameter
#' @param alpha the two-sided confidence interval is level (1-alpha)
#' @param ci_width width of the confidence interval
#'
#' @returns
#' @export
#'
#' @examples
#' pi_cntrl <- 0.03
#' theta <- 0.05
#' rho <- 0.00025
#' out <- design_n_and_r_amplicon_seq(pi_cntrl, theta, rho)
design_n_and_r_amplicon_seq <- function(pi_cntrl, theta, rho, ci_width = 0.01, alpha = 0.05) {
  c <- 1.1
  n <- ceiling(((1/rho) - 1)/(c - 1))
  pi_trt <- pi_cntrl + theta - pi_cntrl * theta
  A_squared <- 1/(1 - pi_cntrl)^4 * ( (1/n + (1 + 1/n) * rho) * ((1 - pi_cntrl)^2 * pi_trt * (1 - pi_trt) + (1 - pi_trt)^2 * pi_cntrl * (1 - pi_cntrl)))
  mult_factor <- qnorm(p = 1 - alpha/2)
  r <- ceiling((2 * mult_factor)^2 * A_squared/ci_width^2)
  # r <- (2 * mult_factor)^2 * A_squared/ci_width^2
  c(n = n, r = r)
}

#' Design n as a function of r in an amplicon-seq experiment
#'
#' Given key parameters of an amplicon-seq study, compute the minimum
#' per-replicate read depth required to attain a target confidence interval
#' width for each candidate replicate count.
#'
#' When `side = two_sided`, `ci_width` is the width of the two-sided interval (i.e., upper limit minus lower limit); when `side = one_sided`, `ci_width` is the upper limit minus the estimate.
#'
#' @param pi_cntrl mutation rate in control condition
#' @param theta editing rate
#' @param rho overdispersion parameter
#' @param ci_width total interval width for `two_sided`, or the margin from
#'   the estimate to the upper confidence bound for `one_sided`
#' @param alpha the confidence level is (1-alpha)
#' @param r_grid integer vector of replicate counts per condition to evaluate
#' @param side `two_sided` (default) or `one_sided`
#'
#' @returns a data frame with one row per value of `r_grid`
#' @export
#'
#' @examples
#' pi_cntrl <- 0.01
#' theta <- 0.0
#' rho <- 0.0005
#' ci_width <- 0.01
#' alpha <- 0.01
#' out <- design_n_and_r_amplicon_seq_v2(pi_cntrl = pi_cntrl, theta = theta,
#' rho = rho, ci_width = ci_width, alpha = alpha)
design_n_and_r_amplicon_seq_v2 <- function(pi_cntrl, theta, rho, ci_width = 0.01, alpha = 0.01, r_grid = seq(1L, 15L), side = c("two_sided", "one_sided")[1]) {
  pi_trt <- pi_cntrl + theta - pi_cntrl * theta
  mult_factor <- if (side == "two_sided") {
    2 * qnorm(1 - alpha / 2)
  } else if (side == "one_sided") {
    qnorm(1 - alpha)
  } else {
    stop("`side` not recognized.")
  }
  B <- 1/(1 - pi_cntrl)^4 * ((1 - pi_cntrl)^2 * pi_trt * (1 - pi_trt) + (1 - pi_trt)^2 * pi_cntrl * (1 - pi_cntrl))
  design_df <- lapply(X = r_grid, FUN = function(r) {
    target <- (ci_width / mult_factor)^2 * r/B
    feasible <- target > rho
    n <- if (feasible) ceiling((1 - rho)/(target - rho)) else Inf
    data.frame(r = r, n = n, feasible = feasible)
  }) |> dplyr::bind_rows()

  return(design_df)
}
