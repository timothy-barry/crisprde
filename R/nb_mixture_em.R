#' Fit a two-component negative binomial mixture by EM
#'
#' Suppose $Y_i$ is drawn from a two-component mixture distribution, i.e.
#' $Y_i = (1 - Z_i) A_i + Z_i B_i,$ where $Z_i \sim \textrm{Bern}(\pi)$,
#' $A_i \sim \textrm{NB}(\mu_a, \theta_a),$ and $B_i \sim \textrm{NB}(\mu_b, \theta_b).$
#' This function estimates the parameters $\pi, \mu_a, \theta_a, \mu_b, \theta_b$ via an EM algorithm.
#' The components are returned in increasing mean order, i.e., the component with the smaller mean is listed first,
#' and that with the larger mean is listed second.
#'
#' Starting means are `mu_multipler_start * mean(y)`; `theta_start` and `pi_start` are used directly. It performs only
#' one optimization, i.e., starting with these initial guesses, it iterates until convergence.
#'
#' Standard M steps maximize the weighted NB likelihood, with size constrained to `[0.001, 100]`.
#' Robust M steps use `fit_rob_nb_univariate()` and its default dispersion bounds, with
#' frequency-weighted memberships and the current component parameters as pilots.
#' Memberships are refreshed after each M step;
#' the mixing proportion retains its usual EM update. Robust updates need not increase the likelihood.
#'
#' An internal optimization is that repeated counts are collapsed and weighted by their frequencies.
#'
#' @param y Nonnegative integer counts, including zeros.
#' @param mu_multipler_start Multipliers of mean(y) for the initial component means.
#' @param theta_start Initial component sizes (variance = mu + mu^2 / theta).
#' @param pi_start Initial probability of component 2.
#' @param maxit Maximum number of EM iterations.
#' @param tol Relative log-likelihood convergence tolerance, or scaled parameter-change
#'   tolerance when `robust_estimation = TRUE`.
#' @param c Tukey tuning constant for both robust mean and dispersion updates.
#' @param robust_estimation Use robust NB fits instead of weighted MLEs in the M step.
#' @return A list containing `converged`, `parameters`, `log_likelihood`, `counts`, `freq`, and
#'   `membership_probabilities`. The named parameter vector contains `mu_lower`,
#'   `mu_upper`, `theta_lower`, `theta_upper`, and `pi`, the higher-mean component's weight.
#'   `counts` contains unique input counts, `freq` their frequencies, and
#'   `membership_probabilities` their higher-mean component membership probabilities.
#'   `log_likelihood` is the ordinary mixture log likelihood, including for robust fits.
#' @examples
#' pi <- 0.3
#' mu_1 <- 5
#' mu_2 <- 10
#' theta_1 <- 0.9
#' theta_2 <- 0.2
#' n <- 5000
#' y_1 <- MASS::rnegbin(n = n, mu = mu_1, theta = theta_1)
#' y_2 <- MASS::rnegbin(n = n, mu = mu_2, theta = theta_2)
#' z <- rbinom(n = n, prob = pi, size = 1)
#' y <- (1 - z) * y_1 + z * y_2
#' fit <- fit_nb_mixture_em_low_level(y)
#' fit <- fit_nb_mixture(y)
fit_nb_mixture_em_low_level <- function(y, mu_multipler_start = c(0.5, 2), theta_start = c(0.5, 0.5), pi_start = 1/3, maxit = 5000L, tol = 1e-6, c = 10, robust_estimation = FALSE) {
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

  mu <- mu_multipler_start * mean(y)
  theta <- theta_start
  pi <- pi_start
  current <- e_step(mu, theta, pi)
  converged <- FALSE
  for (iteration in seq_len(maxit)) {
    previous_loglik <- current$loglik
    previous_parameters <- c(mu, theta, pi)
    weights <- freq * current$tau
    component_n <- colSums(weights)
    pi <- component_n[2] / sum(freq)
    if (robust_estimation) {
      component_fits <- vapply(seq_len(2), function(k) {
        fit_rob_nb_univariate(y = counts, pilot = c(mu = mu[[k]], theta = theta[[k]]),
          weights = weights[, k], c.tukey.beta = c, c.tukey.sigma = c)
      }, c(mu = 0, theta = 0))
      mu <- unname(component_fits["mu", ])
      theta <- unname(component_fits["theta", ])
    } else {
      mu <- colSums(weights * counts) / component_n
      theta <- vapply(seq_len(2), function(k) {
        fit <- stats::optimize(function(log_theta) {
          sum(weights[, k] * stats::dnbinom(counts, mu = mu[k],
            size = exp(log_theta), log = TRUE)) / component_n[k]
        }, interval = c(log(0.001), log(100)), maximum = TRUE, tol = 1e-8)
        exp(fit$maximum)
      }, numeric(1))
    }
    current <- e_step(mu, theta, pi)
    converged <- if (robust_estimation) {
      max(abs(c(mu, theta, pi) - previous_parameters) / (1 + abs(previous_parameters))) <= tol
    } else {
      abs(current$loglik - previous_loglik) <= tol * (1 + abs(previous_loglik))
    }
    if (converged) break
  }

  component_order <- order(mu)
  l <- list(converged = converged,
    parameters = c(mu_lower = mu[component_order[1]], mu_upper = mu[component_order[2]],
      theta_lower = theta[component_order[1]], theta_upper = theta[component_order[2]],
      pi = c(1 - pi, pi)[component_order[2]]),
    log_likelihood = current$loglik,
    membership_probabilities = current$tau[, component_order[2]],
    counts = counts, freq = freq)
  return(l)
}


fit_nb_mixture <- function(y, mu_lower_multiplier_range = c(0.1, 1),
                           mu_upper_multiplier_range = c(1, 10), theta_range = c(0.01, 5),
                           pi_range = c(0.1, 0.9), n_initial_starts = 10,
                           maxit = 5000L, tol = 1e-10, c = 10, robust_estimation = FALSE) {
  # Run standard EM over multiple starting parameters.
  candidate_em_fits <- lapply(X = seq_len(n_initial_starts), FUN = function(i) {
    print(paste0("Running EM algorithm on starting parameter ", i, "."))
    mu_lower_multiplier_start <- runif(n = 1, min = mu_lower_multiplier_range[1], max = mu_lower_multiplier_range[2])
    mu_upper_multiplier_start <- runif(n = 1, min = mu_upper_multiplier_range[1], max = mu_upper_multiplier_range[2])
    mu_multiplier_start <- c(mu_lower_multiplier_start, mu_upper_multiplier_start)
    theta_start <- exp(runif(n = 2, min = log(theta_range[1]), max = log(theta_range[2])))
    pi_start <- runif(n = 1, min = pi_range[1], max = pi_range[2])
    fit <- fit_nb_mixture_em_low_level(y = y,
                                       mu_multipler_start = mu_multiplier_start,
                                       theta_start = theta_start, pi_start = pi_start,
                                       maxit = maxit, tol = tol, c = c,
                                       robust_estimation = FALSE)
  })

  # check for convergence
  converged <- sapply(candidate_em_fits, FUN = function(fit) fit$converged)
  if (any(converged)) {
    retained_fits <- candidate_em_fits[converged]
  } else {
    warning("No EM algorithm fit converged.")
    retained_fits <- candidate_em_fits
  }
  log_likelihoods <- sapply(retained_fits, FUN = function(fit) fit$log_likelihood)
  fit <- retained_fits[[which.max(log_likelihoods)]]
  if (robust_estimation) {
    print("Running robust EM algorithm.")
    parameters <- fit$parameters
    fit <- fit_nb_mixture_em_low_level(y = y,
      mu_multipler_start = unname(parameters[c("mu_lower", "mu_upper")]) / mean(y),
      theta_start = unname(parameters[c("theta_lower", "theta_upper")]), pi_start = parameters[["pi"]],
      maxit = maxit, tol = tol, c = c, robust_estimation = TRUE)
  }
  return(fit)
}
