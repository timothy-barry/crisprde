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
#' Generate synthetic, multivariate amplicon-seq data using an example allele table as input.
#'
#' @param allele_feature_table allele
#' @param frac_mutated fraction of reads containing a mutation
#' @param rho overdispersion parameter
#' @param n_reads number of reads (fixed across replicates)
#' @param beta_mut_length model parameter controlling the relationship between mutation length and allele frequency
#'
#' @returns a count table containing the data
#' @export
#'
#' @examples
#' allele_table_fp <- "/Users/timbarry/research_offsite/external/bauer-lab/rhampseq_bcl11a/crispresso_output_1617/GUIDE0724-africa-edited/CRISPResso_on_1617_OT_0000/Alleles_frequency_table_around_sgRNA_CTAACAGTTGCTTTTATCAC.txt"
#' allele_feature_table <- process_crispresso_allele_table_cas9(allele_table_fp)
#' count_tab <- generate_synthetic_multivariate_amplicon_seq_data(allele_feature_table, n_reads = 100000L)
generate_synthetic_multivariate_amplicon_seq_data <- function(allele_feature_table, use_observed_n_reads = FALSE, sample_size_mu = 100000L, sample_size_theta = 15L, frac_mutated = 0.8, n_rep = 3L, rho = 2e-5, beta_mut_length = -0.1) {
  mod_tab <- allele_feature_table |> dplyr::filter(modified)
  mutation_length <- mod_tab$mutation_length

  # mutation model
  exp_f <- exp(mutation_length * beta_mut_length)
  phi <- exp_f/sum(exp_f)
  tilde_hat <- frac_mutated * phi

  # number of reads
  if (use_observed_n_reads) {
    n_reads <- sum(allele_feature_table$read_count)
  } else {
    n_reads <- MASS::rnegbin(n = n_rep, mu = sample_size_mu, theta = sample_size_theta)
  }

  # simulate reads
  count_tab <- dirmult::simPop(J = n_rep, K = nrow(allele_feature_table), n = n_reads,
                               pi = c(1 - frac_mutated, tilde_hat), theta = rho)$data
  colnames(count_tab) <- c("unmutated", paste0("allele_", seq(1L, nrow(mod_tab))))
  return(count_tab)
}
