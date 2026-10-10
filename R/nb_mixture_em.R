#' Fit a two-component negative binomial mixture by EM
#'
#' Repeated counts are collapsed and weighted by their frequencies. For shifted
#' positive counts, pass `y = umi_count - 1`; returned means are on this scale.
#' Each call uses one initialization. Components are returned in increasing
#' mean order. NB sizes are optimized over log(size) in [-12, 12].
#'
#' @param y Nonnegative integer counts, including zeros.
#' @param mu Initial component means.
#' @param theta Initial component sizes (variance = mu + mu^2 / theta).
#' @param pi Initial probability of component 2.
#' @param maxit Maximum number of EM iterations.
#' @param tol Relative log-likelihood convergence tolerance.
#' @return Component parameters, compressed counts and responsibilities,
#'   log-likelihood and its history, iteration count, and convergence flag.
#' @noRd
fit_nb_mixture_em <- function(y, mu = mean(y) * c(0.5, 2),
                              theta = c(1, 1), pi = 1 / 3,
                              maxit = 10000L, tol = 1e-10) {
  tab <- table(y)
  counts <- as.numeric(names(tab))
  freq <- as.numeric(tab)

  e_step <- function(mu, theta, pi) {
    log_mass <- cbind(
      log1p(-pi) + stats::dnbinom(counts, mu = mu[1], size = theta[1], log = TRUE),
      log(pi) + stats::dnbinom(counts, mu = mu[2], size = theta[2], log = TRUE))
    max_log_mass <- pmax(log_mass[, 1], log_mass[, 2])
    log_mixture <- max_log_mass + log(rowSums(exp(log_mass - max_log_mass)))
    list(tau = exp(log_mass - log_mixture), loglik = sum(freq * log_mixture))
  }

  current <- e_step(mu, theta, pi)
  loglik <- numeric(maxit + 1L)
  loglik[1] <- current$loglik
  converged <- FALSE
  for (iteration in seq_len(maxit)) {
    weights <- freq * current$tau
    component_n <- colSums(weights)
    pi <- component_n[2] / sum(freq)
    mu <- colSums(weights * counts) / component_n
    theta <- vapply(seq_len(2), function(k) {
      fit <- stats::optimize(function(log_theta) {
        sum(weights[, k] * stats::dnbinom(counts, mu = mu[k],
          size = exp(log_theta), log = TRUE)) / component_n[k]
      }, interval = c(-12, 12), maximum = TRUE, tol = 1e-8)
      exp(fit$maximum)
    }, numeric(1))
    current <- e_step(mu, theta, pi)
    loglik[iteration + 1L] <- current$loglik
    converged <- abs(loglik[iteration + 1L] - loglik[iteration]) <=
      tol * (1 + abs(loglik[iteration]))
    if (converged) break
  }

  component_order <- order(mu)
  list(parameters = data.frame(mu = mu[component_order],
      theta = theta[component_order], weight = c(1 - pi, pi)[component_order]),
    counts = data.frame(y = counts, frequency = freq,
      tau_1 = current$tau[, component_order[1]],
      tau_2 = current$tau[, component_order[2]]),
    log_likelihood = current$loglik,
    log_likelihood_trace = loglik[seq_len(iteration + 1L)],
    iterations = iteration, converged = converged)
}
