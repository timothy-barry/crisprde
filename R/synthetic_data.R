#' Generate synthetic amplicon-seq data
#'
#' @param p number of amplicons
#' @param r number of treatment/control samples
#' @param pi_cntrl ([0,1] scalar) background editing rate to use across all control samples
#' @param editing_rate ([0,1] scalar) editing rate in samples where there is nonzero editing
#' @param n_amplicons_nonzero_editing (integer) number of amplicons
#' @param sample_size_mu (positive scalar) mean sample size
#' @param sample_size_theta (positive scalar) size (i.e., overdispersion) parameter for sample size
#' @param beta_binom_rho (positive scalar or vector) dispersion parameter for the beta-binomial distribution
#'
#' @returns
#' @export
#'
#' @examples
#' allele_table_fp <- "/Users/timbarry/research_offsite/external/bauer-lab/rhampseq_bcl11a/crispresso_output_1617/GUIDE0724-africa-edited/CRISPResso_on_1617_OT_0000/Alleles_frequency_table_around_sgRNA_CTAACAGTTGCTTTTATCAC.txt"
#' allele_feature_table <- process_crispresso_allele_table_cas9(allele_table_fp)
generate_synthetic_amplicon_seq_data <- function(p, r, pi_cntrl, editing_rate, n_amplicons_nonzero_editing,
                                                 beta_binom_rho, amplicon_ids, sample_size_mu = 100000L, sample_size_theta = 15L) {
  pi_trt_under_editing <- (pi_cntrl + editing_rate - pi_cntrl * editing_rate)
  pi_cntrl_vect <- rep(pi_cntrl, p)
  pi_trt_vect <- c(rep(pi_trt_under_editing, n_amplicons_nonzero_editing), rep(pi_cntrl, p - n_amplicons_nonzero_editing))
  df_list <- list()
  counter <- 1L
  n_trt_samples <- p * r

  # generate n_mats
  n_mat_list <- lapply(X = seq_len(2L), FUN = function(i) {
    matrix(data = MASS::rnegbin(n = n_trt_samples,
                                mu = sample_size_mu,
                                theta = sample_size_theta), nrow = r, ncol = p)
  })
  n_mat_trt <- n_mat_list[[1]]; n_mat_cntrl <- n_mat_list[[2]]

  # generate k_mats
  if (length(beta_binom_rho) == 1L) beta_binom_rho <- rep(beta_binom_rho, p)
  generate_k_mat <- function(n_mat, pi_vect, beta_binom_rho, p, r) {
    k_mat <- sapply(X = seq_len(p), FUN = function(j) {
      VGAM::rbetabinom(n = r, size = n_mat[,j], prob = pi_vect[j], rho = beta_binom_rho[j])
    })
    return(k_mat)
  }
  k_mat_trt <- generate_k_mat(n_mat_trt, pi_trt_vect, beta_binom_rho, p, r)
  k_mat_cntrl <- generate_k_mat(n_mat_cntrl, pi_cntrl_vect, beta_binom_rho, p, r)

  # set column names
  out_list <- list(amplicon_ids = amplicon_ids,
                   n_mat_trt = n_mat_trt,
                   n_mat_cntrl = n_mat_cntrl,
                   k_mat_trt = k_mat_trt,
                   k_mat_cntrl = k_mat_cntrl,
                   pi_trt_vect = pi_trt_vect,
                   pi_cntrl_vect = pi_cntrl_vect,
                   beta_binom_rho = beta_binom_rho)
  return(out_list)
}


#' Generate synthetic multivariate amplicon-seq data
#'
#' @param n_mutation_types
#' @param n_alleles
#' @param n_covariates
#' @param pi_block
#' @param theta
#' @param phi_block
#' @param gamma
#' @param delta
#'
#' @returns
#' @export
#'
#' @examples
#' n_mutation_types <- 2L
#' n_alleles <- c(50L, 55L)
#' pi_block <- c(0.98, 0.01, 0.01)
#' theta <- 0.75
#' phi_block <- c(0.5, 0.5)
#' gamma_mat <- matrix(c(-1.0, 0.1, -0.8, 0.2), nrow = 2L, byrow = TRUE)
#' delta_mat <- matrix(c(-3, 0.05, -2.75, 0.15), nrow = 2L, byrow = TRUE)
#' rownames(gamma_mat) <- rownames(delta_mat) <- c("mutation_type_1", "mutation_type_2")
#' colnames(gamma_mat) <- colnames(delta_mat) <- c("covariate_1", "covariate_2")
#' covariate_bound <- c(-1, 1)
#' rho <- 0.001
#' n_covariates <- 2L
#' r <- 6L
#' mu_reads <- 50000L
#' theta_reads <- 10
#' sim_data <- generate_synthetic_multivariate_amplicon_seq_data(n_mutation_types = n_mutation_types,
#' n_alleles = n_alleles, n_covariates = n_covariates, pi_block = pi_block,
#' theta = theta, phi_block = phi_block, gamma_mat = gamma_mat, delta_mat = delta_mat,
#' covariate_bound = covariate_bound, rho = rho, r = r, mu_reads = mu_reads,
#' theta_reads = theta_reads)
generate_synthetic_multivariate_amplicon_seq_data <- function(n_mutation_types, n_alleles, n_covariates, pi_block, theta,
                                                              phi_block, gamma_mat, delta_mat, covariate_bound, rho, r,
                                                              mu_reads, theta_reads) {
  # generate the covariate matrix for the different mutation types
  X_list <- lapply(X = seq_len(n_mutation_types), FUN = function(curr_mutation_type) {
    X <- replicate(n = n_covariates, expr = {
      runif(n = n_alleles[curr_mutation_type], min = covariate_bound[1], max = covariate_bound[2]) |> sort()
    })
    colnames(X) <- paste0("covariate_", seq_len(n_covariates))
    rownames(X) <- paste0("type_", curr_mutation_type, "_allele_", seq_len(n_alleles[curr_mutation_type]))
    X
  }) |> setNames(paste0("type_", seq_len(n_mutation_types)))

  # generate background spectrum psi and editing spectrum phi for the different mutation types
  spectra_list <- lapply(X = seq_len(n_mutation_types), FUN = function(curr_mutation_type) {
    # iterate over control and editing spectra
    gamma_coefs <- gamma_mat[curr_mutation_type,]
    delta_coefs <- delta_mat[curr_mutation_type,]
    background_spectrum <- softmax(as.numeric(X_list[[curr_mutation_type]] %*% gamma_coefs))
    editing_spectrum <- softmax(as.numeric(X_list[[curr_mutation_type]] %*% delta_coefs))
    list(background_spectrum = background_spectrum, editing_spectrum = editing_spectrum)
  }) |> setNames(paste0("type_", seq_len(n_mutation_types)))

  # compute background mutation rate vector pi
  pi_mod <- sapply(X = seq_len(n_mutation_types), FUN = function(curr_mutation_type) {
    pi_block[curr_mutation_type + 1L] * spectra_list[[curr_mutation_type]]$background_spectrum
  }) |> setNames(paste0("type_", seq_len(n_mutation_types)))
  pi <- c(pi_block[1], unlist(setNames(pi_mod, NULL)))

  # compute treated mutation rate vector tau
  tau_mod <- sapply(X = seq_len(n_mutation_types), FUN = function(curr_mutation_type) {
    pi_mod[[curr_mutation_type]] + pi[1] * theta * phi_block[curr_mutation_type] * spectra_list[[curr_mutation_type]]$editing_spectrum
  }) |> setNames(paste0("type_", seq_len(n_mutation_types)))
  tau <- c(pi[1] * (1 - theta), unlist(setNames(tau_mod, NULL)))

  # calculate the allele-specific editing rate (for the output)
  theta_tilde <- sapply(X = seq_len(n_mutation_types), FUN = function(curr_mutation_type) {
    theta * phi_block[curr_mutation_type] * spectra_list[[curr_mutation_type]]$editing_spectrum
  }) |> unlist()

  # control read count
  sizes <- MASS::rnegbin(n = r, mu = mu_reads, theta = theta_reads)
  alpha_pi <- get_alpha_from_mu_rho(pi, rho)
  Y_cntrl <- extraDistr::rdirmnom(n = r, size = sizes, alpha = alpha_pi)

  # treated read count
  sizes <- MASS::rnegbin(n = r, mu = mu_reads, theta = theta_reads)
  alpha_tau <- get_alpha_from_mu_rho(tau, rho)
  Y_trt <- extraDistr::rdirmnom(n = r, size = sizes, alpha = alpha_tau)

  # construct the outputs
  Y <- rbind(Y_cntrl, Y_trt)
  X <- do.call(what = rbind, args = X_list)
  allele_df <- data.frame(allele_id = rownames(X),
                          mutation_type = rep(paste0("type_", seq_len(n_mutation_types)),
                                              n_alleles), X)
  allele_df <- dplyr::bind_rows(data.frame(allele_id = "unmutated",
                                           mutation_type = NA_character_), allele_df)
  rownames(allele_df) <- NULL
  rep_id <- paste0("rep_", seq_len(2 * r))
  covariate_df <- data.frame(rep_id = rep_id,
                             treated = c(rep(FALSE, r), rep(TRUE, r)))
  rownames(Y) <- rep_id
  colnames(Y) <- c("unmutated", allele_df$allele_id)
  l <- list(Y = Y, X = X, allele_df = allele_df, covariate_df = covariate_df, theta_tilde = theta_tilde, pi = pi, tau = tau)
  return(l)
}


# softmax helper function
softmax <- function(x) exp(x)/sum(exp(x))
# DM helper function -- get alpha from mu, rho
get_alpha_from_mu_rho <- function(mu, rho) mu * (1 - rho) / rho
