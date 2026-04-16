#' Compute basic statistics on multivariate amplicon-seq data
#'
#' This function computes the editing rate and attributable fraction of each allele in a multivariate allele table.
#'
#' @param count_matrix the allele count matrix
#' @param covariate_df the allele covariate data frame
#'
#' @returns a data frame with columns `editing_rate` and `attributable_fraction`
#' @export
#'
#' @examples
#' data_list <- readRDS("/Users/timbarry/research_offsite/external/bauer-lab/rhampseq_bcl11a/multivariate_count_tables/combined/CRISPResso_on_1450_OT_0000.rds")
#' count_matrix <- data_list$count_matrix
#' covariate_df <- data_list$covariate_df
compute_basic_statistics_on_multivariate_amplicon_seq_data <- function(count_matrix, covariate_df) {
  # compute theta_tilde and attributable fraction
  rownames(covariate_df) <- covariate_df$replicate_id
  treated <- covariate_df[rownames(count_matrix), "treated"]
  cntrl_mat <- count_matrix[!treated,]
  treated_mat <- count_matrix[treated,]
  mutated_alleles <- colnames(cntrl_mat)[colnames(cntrl_mat) != "unmutated"]

  p_0_cntrl <- sum(cntrl_mat[,"unmutated"])/sum(cntrl_mat)
  p_1_trt <- Matrix::colSums(treated_mat[,mutated_alleles])/sum(treated_mat)
  p_1_cntrl <- Matrix::colSums(cntrl_mat[,mutated_alleles])/sum(cntrl_mat)
  theta_tilde <- pmax((p_1_trt - p_1_cntrl)/p_0_cntrl, 0)
  attributable_fraction <- pmax(ifelse(p_1_trt < 1e-12, NA_real_, (p_1_trt - p_1_cntrl)/p_1_trt), 0)

  out_df <- data.frame(allele_id = mutated_alleles,
                       theta_tilde = theta_tilde,
                       attributable_fraction = attributable_fraction)
  rownames(out_df) <- NULL
  out_df |> dplyr::arrange(-theta_tilde)
}


#' Collapse allele table
#'
#' @param count_matrix an allele count matrix
#' @param allele_df the allele data frame
#' @param bucket_breaks integer vector indicating the levels to break the mutation lengths into
#' @param bucket_labels character vector indicating the bucket labels
#'
#' @returns the collapsed matrix
#' @export
#'
#' @examples
#' data_list <- readRDS("/Users/timbarry/research_offsite/external/bauer-lab/rhampseq_bcl11a/multivariate_count_tables/combined/CRISPResso_on_1450_OT_0000.rds")
#' count_matrix <- data_list$count_matrix
#' allele_df <- data_list$allele_df
#' covariate_df <- data_list$covariate_df
#' collapsed_matrix_list <- collapse_count_matrix_by_allele(count_matrix, allele_df)
#' collapsed_res <- compute_basic_statistics_on_multivariate_amplicon_seq_data(collapsed_matrix_list$collapsed_matrix, covariate_df)
#'
collapse_count_matrix_by_allele <- function(count_matrix, allele_df,
                                            bucket_breaks = c(seq(0L, 49L), Inf),
                                            bucket_labels = c(seq(1L, 49L), "50+")) {
  unmutated_count_vect <- as.matrix(count_matrix[, "unmutated",drop=FALSE])
  replicate_levels <- rownames(count_matrix)
  allele_df_sub <- allele_df |> dplyr::filter(allele_id != "unmutated")
  allele_df_sub$bucketed_length <- cut(x = allele_df_sub$mutation_length, breaks = bucket_breaks, labels = bucket_labels)

  # get sparse triplet matrix representation
  count_matrix_t <- as(count_matrix[, colnames(count_matrix) != "unmutated"], Class = "TsparseMatrix")
  count_matrix_df <- data.frame(i = count_matrix_t@i, j = count_matrix_t@j, x = count_matrix_t@x)
  count_matrix_df$allele_id  <- factor(colnames(count_matrix_t)[count_matrix_df$j + 1L])
  count_matrix_df$replicate_id <- factor(rownames(count_matrix_t)[count_matrix_df$i + 1L])
  count_matrix_df <- dplyr::left_join(x = count_matrix_df, y = allele_df_sub, by = c("allele_id"))

  # combine counts across mutation buckets via dplyr
  count_matrix_df_sum <- count_matrix_df |>
    dplyr::group_by(replicate_id, mutation_type, bucketed_length) |>
    dplyr::summarize(read_count = sum(x)) |>
    dplyr::ungroup() |>
    dplyr::arrange(mutation_type, bucketed_length) |>
    dplyr::mutate(replicate_id = as.character(replicate_id), mutation_bucket = paste0(mutation_type, "_", bucketed_length))

  # reconstruct sparse matrix
  mutation_bucket_levels <- unique(count_matrix_df_sum$mutation_bucket)
  count_matrix_df_sum$mutation_bucket_idx <- match(count_matrix_df_sum$mutation_bucket, mutation_bucket_levels)
  count_matrix_df_sum$replicate_idx <- match(count_matrix_df_sum$replicate_id, replicate_levels)
  count_matrix <- Matrix::sparseMatrix(
    i = count_matrix_df_sum$replicate_idx,
    j = count_matrix_df_sum$mutation_bucket_idx,
    x = count_matrix_df_sum$read_count,
    dims = c(length(replicate_levels), length(mutation_bucket_levels))
  ) |> as.matrix()
  rownames(count_matrix) <- replicate_levels
  colnames(count_matrix) <- mutation_bucket_levels

  # append the unmutated column
  collapsed_matrix <- cbind(unmutated_count_vect, count_matrix)

  # construct the collapsed allele data frame
  mutation_bucket_split <- strsplit(x = colnames(count_matrix), split = "_", fixed = TRUE)
  mutation_type <- sapply(mutation_bucket_split, FUN = function(l) l[[1]])
  mutation_length <- sapply(mutation_bucket_split, FUN = function(l) l[[2]])
  mutation_length_int <- gsub(pattern = "+", replacement = "", x = mutation_length, fixed = TRUE) |> as.integer()
  collapsed_allele_df <- data.frame(allele_id = mutation_bucket_levels,
                                    mutation_type = mutation_type,
                                    mutation_length = mutation_length_int)
  collapsed_allele_df <- rbind(data.frame(allele_id = "unmutated", mutation_type = NA_character_, mutation_length = NA_integer_),
                               collapsed_allele_df)

  # construct output
  ret <- list(collapsed_matrix = collapsed_matrix, collapsed_allele_df = collapsed_allele_df)

  # collapse the allele feature table
  return(ret)
}


#' Fit multivariate Bayesian regression model
#'
#' @returns
#' @export
#'
#' @examples
#' # prepare data input
#' data_list <- readRDS("/Users/timbarry/research_offsite/external/bauer-lab/rhampseq_bcl11a/multivariate_count_tables/combined/CRISPResso_on_1450_OT_0000.rds")
#' count_matrix <- data_list$count_matrix
#' allele_df <- data_list$allele_df
#' covariate_df <- data_list$covariate_df
#' collapsed_matrix_list <- collapse_count_matrix_by_allele(count_matrix, allele_df)
#' count_matrix <- collapsed_matrix_list$collapsed_matrix
#' allele_df <- collapsed_matrix_list$collapsed_allele_df
#' # for now -- restrict attention to deletions
#' to_keep <- allele_df$mutation_type %in% c(NA, "deletion")
#' count_matrix <- count_matrix[,to_keep]
#' allele_df <- allele_df[to_keep,]
#'
#' # get rho from univariate analysis
#' dat_list <- readRDS("/Users/timbarry/research_offsite/external/bauer-lab/rhampseq_bcl11a/marginal_count_tables/data_list_by_grna.rds")[["1450"]]
#' univariate_res <- run_freqentist_amplicon_seq_analysis(data_list = dat_list)
#' rho_hat <- univariate_res$result_df |> dplyr::filter(amplicon_id == "1450_OT_0000") |> dplyr::pull(rho_hat)
#'
#' # formula object
#' formula_object <- formula(~ mutation_length + 0)
#'
#' # hyperparameters
#' # pi (control mean)
#' mu_pi_block <- c(0.999, 0.001)
#' kappa_pi_block <- 100
#'
#' # theta (editing rate)
#' mu_theta <- 0.5
#' kappa_theta <- 20
#'
#' # block level editing
#' mu_phi_block <- 1
#' kappa_phi_block <- 1
#'
#' # control side regression
#' mu_gamma <- matrix(-1)
#' sigma_gamma <- matrix(0.5)
#'
#' # treatment side regression
#' mu_delta <- matrix(-1)
#' sigma_delta <- matrix(0.5)
#'
#' # fit model
#' fit <- fit_multivariate_bayesian_regression_model(count_matrix, allele_df, covariate_df, rho_hat, formula_object,
#' mu_pi_block, kappa_pi_block,
#' mu_theta, kappa_theta,
#' mu_phi_block, kappa_phi_block,
#' mu_gamma, sigma_gamma,
#' mu_delta, sigma_delta)
fit_multivariate_bayesian_regression_model <- function(count_matrix, allele_df, covariate_df, rho_hat, formula_object,
                                                       mu_pi_block, kappa_pi_block,
                                                       mu_theta, kappa_theta,
                                                       mu_phi_block, kappa_phi_block,
                                                       mu_gamma, sigma_gamma,
                                                       mu_delta, sigma_delta) {
  # construct the model matrix
  allele_df <- allele_df |> dplyr::filter(allele_id != "unmutated")
  X <- model.matrix(object = formula_object, data = allele_df) |> scale()

  # set key model parameters
  unique_mutation_types <- unique(allele_df$mutation_type)
  m <- length(unique_mutation_types)
  start_stop_info <- sapply(X = unique_mutation_types, FUN = function(curr_type) {
    range(which(allele_df$mutation_type == curr_type))
  }) |> t()
  type_start <- start_stop_info[,1]
  type_end <- start_stop_info[,2]
  q_t <- apply(X = start_stop_info, FUN = function(r) r[2] - r[1] + 1L, MARGIN = 1L)
  q <- nrow(X)
  p <- ncol(X)
  Y_trt <- count_matrix[covariate_df$replicate_id,][covariate_df$treated,,drop=FALSE]
  Y_cntrl <- count_matrix[covariate_df$replicate_id,][!covariate_df$treated,,drop=FALSE]
  storage.mode(Y_cntrl) <- "integer"
  storage.mode(Y_trt) <- "integer"
  r_trt <- nrow(Y_trt)
  r_cntrl <- nrow(Y_cntrl)
  rho <- rho_hat

  # organize hyperparameters
  alpha_beta_theta <- beta_mu_kappa_to_alpha_beta(mu = mu_theta, kappa = kappa_theta)
  alpha_theta <- alpha_beta_theta[["alpha"]]
  beta_theta <- alpha_beta_theta[["beta"]]

  # initialize the model
  stan_file <- "~/research_code/crispr_safe_parent/crisprde/inst/stan/reg_model.stan" # system.file("stan", "reg_model.stan", package = "crisprde")
  model <- cmdstanr::cmdstan_model(stan_file)
  stan_data <- list(
    m = m,
    q = q,
    p = p,
    r_cntrl = r_cntrl,
    r_trt = r_trt,
    Y_cntrl = Y_cntrl,
    Y_trt = Y_trt,
    rho = rho_hat,
    q_t = q_t,
    type_start = type_start,
    type_end = type_end,
    X = X,
    mu_pi_block = mu_pi_block,
    kappa_pi_block = kappa_pi_block,
    alpha_theta = alpha_theta,
    beta_theta = beta_theta,
    mu_phi_block = mu_phi_block,
    kappa_phi_block = kappa_phi_block,
    mu_gamma = mu_gamma,
    sigma_gamma = sigma_gamma,
    mu_delta = mu_delta,
    sigma_delta = sigma_delta
  )

  # fit model
  fit <- model$sample(
    data = stan_data,
    seed = 4,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95,
    max_treedepth = 12,
    refresh = 250
  )

  return(fit)
}
