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
