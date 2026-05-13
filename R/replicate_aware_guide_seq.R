#' Truncated Bernoulli product (TBP) distribution
#'
#' Simulate draws from a TBP distribution
#'
#' @param m number of samples to draw
#' @param pi the TBP parameter vector of length r
#'
#' @returns a matrix of dimension r x m of TBP draws
#' @examples
#' m <- 5000L
#' pi <- c(0.1, 0.05, 0.1)
#' X <- r_tbp(m, pi)
r_tbp <- function(m, pi) {
  Omega <- pmf_tbp(pi)
  s <- sample(x = seq_len(nrow(Omega)), size = m, replace = TRUE, prob = Omega$pmf)
  out <- Omega[s,seq_along(pi)] |> as.matrix() |> t()
  colnames(out) <- NULL
  return(out)
}

# helper function to generate the Omega matrix of combinations
generate_omega <- function(r) {
  Omega <- rep(list(c(0L, 1L)), r) |> setNames(paste0("x_", seq_len(r))) |> expand.grid()
  Omega <- Omega[rowSums(Omega) >= 1L,]
  rownames(Omega) <- NULL
  return(Omega)
}

# helper function to convolve two pmfs
convolve_pmfs <- function(a, b) {
  pmf <- convolve(a, rev(b), type = "open")
  return(pmf)
}

# helper function to convolve a list of pmfs
convolve_pmf_list <- function(pmf_list) {
  if (length(pmf_list) == 1L) {
    out <- pmf_list[[1]]
  } else {
    out <- Reduce(convolve_pmfs, pmf_list)
  }
  out <- pmax(out, 1e-50)
  out <- out / sum(out)
  return(out)
}

#' Returns the pmf of a truncated Bernoulli product (TBP) distribution
#'
#' @param pi parameter of the TBP distribution
#'
#' @returns a data frame with columns (x_1, x_2, \dots, x_r, pmf), where pmf gives the probability of a given (x_1, x_2, \dots, x_r) vector
#' @examples
#' Omega <- pmf_tbp(c(0.1, 0.05, 0.02))
#' Omega <- pmf_tbp(c(0.4, 0.6, 0.2))
pmf_tbp <- function(pi) {
  r <- length(pi)
  Omega <- generate_omega(r)
  denom <- 1 - prod(1 - pi)
  pmf <- apply(X = Omega, MARGIN = 1, FUN = function(x) {
    prod(ifelse(x, pi, 1 - pi))
  })/denom
  Omega$pmf <- pmf
  return(Omega)
}


#' Fit truncated Bernoulli product (TBP) model via MLE
#'
#' Fits a TBP model to an r x m matrix of occupancies.
#'
#' Returns NULL if the regularity conditions fail to hold.
#'
#' @param X binary occupancy matrix of dimension r (number of replicates) by m (number of windows)
#'
#' @returns the fitted MLE pi-hat (or NULL if the regularity conditions fail to hold)
#' @examples
#' m <- 5000L
#' pi <- c(0.2, 0.12, 0.04)
#' X <- r_tbp(m, pi)
#' pi_hat <- fit_tbp_model(X)
fit_tbp_model <- function(X) {
  # unconstrained estimate
  q_hat <- rowMeans(X)
  # verify regularity conditions
  if (sum(q_hat) > 1 && all(q_hat > 0) && all(q_hat < 1)) {
    g <- function(c) 1 - prod(1 - c * q_hat) - c
    root_res <- uniroot(g, lower = 1e-10, upper = 1 - 1e-10)
    pi_hat <- root_res$root * q_hat
  } else {
    pi_hat <- NULL
  }
  return(pi_hat)
}


#' shifted negative binomial (SNB) distribution
#'
#' Sample m_plus observations from an SNB distribution with parameters (mu, theta)
#'
#' @param m_plus the number of samples to generate (typically equal to the number of occupied windows)
#' @param mu mean parameter
#' @param theta size parameter
#'
#' @returns a vector of snb variates
#' @examples
#' mu <- 10
#' theta <- 0.5
#' m_plus <- 1000L
#' y_plus <- r_snb(m_plus, mu, theta)
r_snb <- function(m_plus, mu, theta) {
  MASS::rnegbin(n = m_plus, mu = mu, theta = theta) + 1L
}


#' Simulate multi-replicate guide-seq data
#'
#' @param pi parameter vector of TBP model
#' @param mu_vect vector of mu parameters for SNB models
#' @param theta_vect vector of theta parameters for SNB models
#'
#' @returns
#' @examples
#' # NULL DATA
#' pi <- c(0.05, 0.1, 0.02)
#' mu_vect <- c(10, 6, 15)
#' theta_vect <- c(2, 5, 0.3)
#' m <- 10000
#' null_dat <- simulate_multirep_guideseq_data(pi, mu_vect, theta_vect, m)
#'
#' # ALTERNATIVE DATA
#' pi <- c(0.5, 0.6, 0.4)
#' mu_vect <- c(80, 200, 50)
#' theta_vect <- c(20, 21, 15)
#' m <- 15
#' alt_dat <- simulate_multirep_guideseq_data(pi, mu_vect, theta_vect, m)
#'
#' # COMBINED DATA
#' Y_mat <- cbind(alt_dat, null_dat)
#' colnames(Y_mat) <- paste0("window_", seq_len(ncol(Y_mat)))
simulate_multirep_guideseq_data <- function(pi, mu_vect, theta_vect, m) {
  X <- r_tbp(m, pi)
  Y <- sapply(X = seq_along(pi), FUN = function(i) {
    y_plus <- r_snb(m_plus = sum(X[i,]), mu = mu_vect[i], theta = theta_vect[i])
    y <- integer(m)
    y[X[i,] == 1L] <- y_plus
    return(y)
  }) |> t()
  return(Y)
}


#' Run multivariate guide-seq method
#'
#' Runs the multivariate guide-seq method on a replicate-by-window matrix of counts
#'
#' @param Y_mat the r x m integer matrix of UMI counts, where r is the number of replicates and m is the number of windows
#' @param incorporate_occupancy_info a boolean (T/F) indicating whether to incorporate occupancy information into the p-value calculation
#'
#' @returns a data frame containing a p-value for each window
#' @examples
#' set.seed(42)
#' # NULL DATA
#' pi <- c(0.05, 0.1, 0.02)
#' mu_vect <- c(10, 6, 15)
#' theta_vect <- c(2, 5, 0.3)
#' m <- 10000
#' null_dat <- simulate_multirep_guideseq_data(pi, mu_vect, theta_vect, m)
#'
#' # ALTERNATIVE DATA
#' pi <- c(0.5, 0.6, 0.4)
#' mu_vect <- c(80, 200, 50)
#' theta_vect <- c(20, 21, 15)
#' m_alt <- 15
#' alt_dat <- simulate_multirep_guideseq_data(pi, mu_vect, theta_vect, m_alt)
#'
#' # COMBINED DATA
#' Y_mat <- cbind(alt_dat, null_dat)
#' colnames(Y_mat) <- paste0("window_", seq_len(ncol(Y_mat)))
#' incorporate_occupancy_info <- TRUE
#' multiplicity_alpha <- 0.2
#'
#' # RUN METHOD
#' res_df <- run_multireplicate_guideseq_method(Y_mat, incorporate_occupancy_info = TRUE)
#'
#' # EVALUATE RESULT
#' n_correct_nominations <- sum(res_df$nominated_window[seq(1, m_alt)])
#' n_total_nominations <- sum(res_df$nominated_window)
#' fdp <- (n_total_nominations - n_correct_nominations)/n_total_nominations
run_multireplicate_guideseq_method <- function(Y_mat, incorporate_occupancy_info = TRUE,
                                               multiplicity_alpha = 0.2, lambda = NULL, c_tukey_beta = 5,
                                               c_tukey_sigma = 5, robust_fit = TRUE) {
  if (is.null(colnames(Y_mat))) warning("Y_mat must have column names (to identify the windows).")
  X <- Y_mat > 0
  storage.mode(X) <- "integer"
  Omega <- as.matrix(generate_omega(nrow(Y_mat)))

  # 1. fit occupancy model and compute marginal occupancy probabilities
  if (incorporate_occupancy_info) {
    pi_hat <- fit_tbp_model(X)
    if (!is.null(pi_hat)) {
      # get pmf of fitted tbp distribution
      tbp_pattern_df <- pmf_tbp(pi_hat)
      if (is.null(lambda)) {
        lambda <- -median(colSums(Y_mat) - colSums(X))/(sum(log(pi_hat)))
      }
    } else {
      warning("Cannot fit occupancy model; defaulting to count-only model.")
      incorporate_occupancy_info <- FALSE
    }
  }
  col_keys <- apply(X, 2, paste0, collapse = "")
  if (!incorporate_occupancy_info) {
    omega_keys <- apply(Omega, 1, paste0, collapse = "")
    occupancy_pattern_map <- match(col_keys, omega_keys)
    lambda <- NA
  }

  # 2. compute shifted NB models and compute convolved pms and cdfs over all combinations
  mu_theta_hat_mat <- apply(X = Y_mat, MARGIN = 1, FUN = function(curr_row) {
    y_plus <- curr_row[curr_row > 0]
    if (robust_fit) {
      fit <- fit_rob_nb_univariate(y = y_plus - 1, c.tukey.beta = c_tukey_beta, c.tukey.sigma = c_tukey_sigma)
    } else {
      fit <- fit_nb_univariate(y = y_plus - 1)
    }
    fit[c("mu", "theta")]
  }) |> t()
  # generate pmfs
  pmf_list <- apply(X = mu_theta_hat_mat, MARGIN = 1, FUN = function(curr_row) {
    max_count <- qnbinom(p = 1e-50, mu = curr_row[["mu"]], size = curr_row[["theta"]], lower.tail = FALSE)
    dnbinom(x = seq(0L, max_count), mu = curr_row[["mu"]], size = curr_row[["theta"]])
  }, simplify = FALSE)
  # compute all convolution combinations
  conv_pmf_list <- apply(X = Omega, MARGIN = 1, FUN = function(curr_row) {
    convolve_pmf_list(pmf_list[as.logical(curr_row)])
  }, simplify = FALSE)
  # compute right-tail probability list
  right_tail_prob_list <- lapply(X = conv_pmf_list, FUN = function(curr_pmf) {
    rev(cumsum(rev(curr_pmf)))
  })

  # 3. compute the p-value for each window
  p_val_and_test_stat_mat <- sapply(X = seq_len(ncol(Y_mat)), FUN = function(j) {
    y <- Y_mat[,j]
    total_umi_count <- sum(y)
    if (!incorporate_occupancy_info) {
      occupancy_pattern_idx <- occupancy_pattern_map[j]
      test_stat <- total_umi_count - sum(X[,j])
      right_tail_prob_vector <- right_tail_prob_list[[occupancy_pattern_idx]]
      p_val <- get_p_value_given_test_stat_prob_vector(test_stat, right_tail_prob_vector)
    } else {
      # compute test stat
      test_stat <- (total_umi_count - sum(X[,j])) - lambda * sum(X[,j] * log(pi_hat))
      p_val <- sapply(X = seq_len(nrow(tbp_pattern_df)), FUN = function(i) {
        curr_pmf <- right_tail_prob_list[[i]]
        sum_start <- ceiling(test_stat + lambda * sum(Omega[i,] * log(pi_hat)))
        nb_piece <- get_p_value_given_test_stat_prob_vector(test_stat = sum_start,
                                                            right_tail_prob_vector = curr_pmf)
        tbp_piece <- tbp_pattern_df$pmf[i]
        return(tbp_piece * nb_piece)
      }) |> sum()
    }
    return(c(test_stat, p_val))
  }) |> t()

  # 4. multiplicity correction
  test_stats <- p_val_and_test_stat_mat[,1]
  p_vals <- pmin(1, p_val_and_test_stat_mat[,2])
  q_vals <- p.adjust(p = p_vals, method = "BH")
  umi_counts <- colSums(Y_mat)
  nominated_window <- (q_vals < multiplicity_alpha)

  res_df <- data.frame(window = colnames(Y_mat),
                       p_value = p_vals,
                       test_stat = test_stats,
                       nominated_window = nominated_window,
                       umi_count = umi_counts,
                       lambda = lambda,
                       occupancy_pattern = col_keys) |> dplyr::arrange(p_value)
  rownames(res_df) <- NULL
  ests_list <- list(mu_theta_hat_mat = mu_theta_hat_mat)
  if (incorporate_occupancy_info) ests_list$pi_hat <- pi_hat
  ret <- list(res_df = res_df, ests_list = ests_list)
  return(ret)
}


get_p_value_given_test_stat_prob_vector <- function(test_stat, right_tail_prob_vector) {
  if (test_stat < 0) {
    p_val <- 1
  } else {
    idx <- test_stat + 1L
    if (idx > length(right_tail_prob_vector)) {
      p_val <- right_tail_prob_vector[length(right_tail_prob_vector)]
    } else {
      p_val <- right_tail_prob_vector[idx]
    }
  }
  return(p_val)
}


#' Cluster loci
#'
#' Clusters loci via single-linkage clustering
#'
#' @param count_df a data frame containing columns chr, coord
#' @param thresh clustering threshold; occupied bases within this distance are clustered together
#' @param padding amount of padding to add to either side of each cluster
#'
#' @returns `count_df` with the additional columns appended:
#' - `group_string` (a string indicating the cluster to which a given base belongs)
#' - `group_chr` (coordinate of the group)
#' - `min_group_coord` (minumum coordinate of the group)
#' - `max_group_coord` (maximum coordinate of the group)
#' @export
#'
#' @examples
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' count_df <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#'  dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & treated & replicate_id %in% 1:2) |>
#'  dplyr::filter(chr != "chrM") |>
#'  dplyr::select(chr, coord, strand, umi_count, primer_type, replicate_id)
#' clustered_count_df <- cluster_loci(count_df)
cluster_loci <- function(count_df, thresh = 100L, padding = 25L) {
  # 1. simple distance
  count_df_w_dist <- count_df |>
    dplyr::group_by(chr) |>
    dplyr::arrange(chr, coord) |>
    dplyr::mutate(simple_dist = c(NA, diff(coord)))
  curr_group_id <- 0L
  ds <- count_df_w_dist$simple_dist
  group_id <- integer(length = nrow(count_df))
  for (i in seq(1, nrow(count_df))) {
    if (ds[i] > thresh || is.na(ds[i])) {
      curr_group_id <- curr_group_id + 1
    }
    group_id[i] <- curr_group_id
  }
  count_df_w_dist_and_group_string <- count_df_w_dist |>
    dplyr::ungroup() |>
    dplyr::mutate(group_id = group_id) |>
    dplyr::group_by(group_id) |>
    dplyr::mutate(group_chr = chr[1],
                  min_group_coord = min(coord) - padding,
                  max_group_coord = max(coord) + padding) |>
    dplyr::ungroup() |>
    dplyr::mutate(group_string = factor(paste0(group_chr, ":", min_group_coord, "-", max_group_coord)),
                  group_id = NULL)
  return(count_df_w_dist_and_group_string)
}


#' Construct replicate count table
#'
#' Takes the output from the GUIDE-seq pipeline as input; outputs a count matrix for statistical modeling
#'
#' @param count_df a data frame with columns `chr`, `coord`, `strand`, `replicate_id`, and `umi_count`.
#' @param thresh clustering threshold; bases within this distance are clustered together into a window
#' @param padding amount of padding to add to either end of each window
#'
#' @returns an integer matrix with replicates in the rows and windows in the columns. An entry of the matrix indicates the number of UMIs observed within a given window and replicate.
#' @export
#'
#' @examples
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' count_df <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#'  dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & treated & replicate_id %in% 1:2, chr != "chrM") |>
#'  dplyr::select(chr, coord, strand, umi_count, primer_type, replicate_id)
#' clustered_count_df <- cluster_loci(count_df)
#' Y_mat <- construct_replicate_count_table(clustered_count_df)
construct_replicate_count_table <- function(clustered_count_df) {
  # if primer_type is present, append that to rep-id
  if ("primer_type" %in% colnames(clustered_count_df)) {
    clustered_count_df <- clustered_count_df |> dplyr::mutate(replicate_id = paste0(replicate_id, "-", primer_type), primer_type = NULL)
  }
  # sum over UMIs within a given (window, replicate) pair
  collapsed_count_df <- clustered_count_df |>
    dplyr::group_by(replicate_id, group_string) |>
    dplyr::summarize(umi_count = sum(umi_count), .groups = "drop") |>
    dplyr::arrange(replicate_id, group_string)

  # construct Y_mat
  replicate_idx <- match(x = collapsed_count_df$replicate_id, unique(collapsed_count_df$replicate_id))
  group_idx <- match(x = collapsed_count_df$group_string, unique(collapsed_count_df$group_string))
  Y_mat <- Matrix::sparseMatrix(i = replicate_idx,
                                j = group_idx,
                                x = collapsed_count_df$umi_count) |>
    as.matrix()
  rownames(Y_mat) <- unique(collapsed_count_df$replicate_id)
  colnames(Y_mat) <- unique(collapsed_count_df$group_string)
  return(Y_mat)
}


#' Tune hyperparameters
#'
#' @param Y_mat_trt treated count matrix
#' @param Y_mat_cntrl control count matrix
#' @param c_grid robust hyperparm grid
#' @param lambda_grid lambda grid
#' @param multiplicity_alpha nominal fdr
#' @param max_false_discs maximum false discoveries permitted in the control condition
#'
#' @returns a list with elements `selected_params`, `selected_trt_run`, `selected_cntrl_run`, `grid_results`
#' @export
#'
#' @examples
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' count_df_all <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#' dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & replicate_id %in% 1:2, chr != "chrM")
#' Y_mat_trt <- count_df_all |> dplyr::filter(treated) |> cluster_loci() |> construct_replicate_count_table()
#' Y_mat_cntrl <- count_df_all |> dplyr::filter(treated) |> cluster_loci() |> construct_replicate_count_table()
#'
#' hyperparam_out <- tune_hyperparameters(Y_mat_trt = Y_mat_trt, Y_mat_cntrl = Y_mat_cntrl, c_grid = c(5, 25), lambda_grid = c(5, 50))
tune_hyperparameters <- function(Y_mat_trt, Y_mat_cntrl, c_grid = c(1, 5, 10, 25, 50, 100),
                                 lambda_grid = c(0, 5, 10, 20, 50, 100, 500), multiplicity_alpha = 0.2, max_false_discs = 1L) {
  condition_grid <- c("trt", "cntrl")
  grid <- expand.grid(c = c_grid, lambda = lambda_grid, condition = condition_grid)

  # helper function to run method on one parameter configuration
  run_one_grid_row <- function(i) {
    curr_row <- grid[i, , drop = FALSE]
    curr_condition <- as.character(curr_row$condition)
    curr_c <- curr_row$c
    curr_lambda <- curr_row$lambda
    curr_Y_mat <- if (curr_condition == "trt") Y_mat_trt else Y_mat_cntrl

    fit_res <- run_multireplicate_guideseq_method(Y_mat = curr_Y_mat, incorporate_occupancy_info = TRUE, robust_fit = TRUE,
                                                  lambda = curr_lambda, c_tukey_beta = curr_c, c_tukey_sigma = curr_c,
                                                  multiplicity_alpha = multiplicity_alpha)
    list(params = curr_row, res = fit_res)
  }
  mc_cores <- max(1L, parallel::detectCores() - 1L)

  # apply method to analyze data
  grid_results <- parallel::mclapply(
    X = seq_len(nrow(grid)),
    FUN = run_one_grid_row,
    mc.cores = mc_cores
  )

  summary_df <- lapply(X = grid_results, FUN = function(curr_res) {
    curr_res$params |>
      dplyr::mutate(n_discoveries = sum(curr_res$res$res_df$nominated_window))
  }) |> data.table::rbindlist() |>
    tidyr::pivot_wider(names_from = "condition", values_from = n_discoveries)
  if (any(summary_df$cntrl <= max_false_discs)) {
    selected_params <- summary_df |>
      dplyr::filter(cntrl <= max_false_discs) |>
      dplyr::arrange(dplyr::desc(trt), lambda, c) |>
      dplyr::slice(1)
    trt_idx <- sapply(grid_results, FUN = function(curr_res) {
      curr_res$params$c == selected_params$c && curr_res$params$lambda == selected_params$lambda && curr_res$params$condition == "trt"
    }) |> which()
    cntrl_idx <- sapply(grid_results, FUN = function(curr_res) {
      curr_res$params$c == selected_params$c && curr_res$params$lambda == selected_params$lambda && curr_res$params$condition == "cntrl"
    }) |> which()
    selected_trt_run <- grid_results[[trt_idx]]$res
    selected_cntrl_run <- grid_results[[cntrl_idx]]$res
  } else {
    selected_params <- NA
    selected_trt_run <- NA
    selected_cntrl_run <- NA
  }
  ret <- list(selected_params = selected_params,
              selected_trt_run = selected_trt_run,
              selected_cntrl_run = selected_cntrl_run,
              grid_results = grid_results, summary_df = summary_df)
  return(ret)
}


#' Load CRISPRitz output
#'
#' @param targets_file_path file path to the targets.txt file outputted by CRISPRitz
#'
#' @returns a data frame containing the CRISPRitz target output
#' @export
#'
#' @examples
#' homology_df <- load_crispritz_output("/Users/timbarry/research_offsite/external/bauer-lab/guideseq_elane/crispritz/crispritz_CCCCGGCAGAAACGTCCGCG.hg38.targets.txt")
load_crispritz_output <- function(targets_file_path) {
  df <- readr::read_delim(file = targets_file_path) |>
    dplyr::rename("bulge_type" = "#Bulge type", "n_mismatches" = "Mismatches",
                  "n_bulges" = "Bulge Size", "n_total_changes" = "Total",
                  "gRNA" = "crRNA", "dna" = "DNA", "chromosome" = "Chromosome",
                  "posit" = "Position", "cluster_posit" = "Cluster Position",
                  "strand" = "Direction")
  dna_width <- nchar(gsub(pattern = "-", replacement = "", x = df$dna)) - 3L # subtract the PAM
  df$dna_width <- dna_width
  return(df)
}


#' Overlap homology and result df
#'
#' Add colummns:
#' - `has_homology_hit` (indicating whether there is a homology hit inside the window)
#' - `homology_n_mismatches` (indicating the number of mismatches between aligned spacer and protospacer sequence)
#' - `homology_n_bulges` (indicating number of bulges between aligned spacer and protospacer sequence)
#' - `homology_posit` (indicating the position of the start of the protospacer)
#' - `homology_strand` (indicating whether the protospacer is on the plus or minus strand)
#' - `homology_dna` (aligned protospacer sequence)
#' - `homology_gRNA` (aligned spacer sequence)
#'
#' @param result_df output of run_multireplicate_guideseq_method
#' @param homology_df a homology data frame (for instance, output of load_crispritz_output)
#'
#' @returns an augmented version of `result_df` with homology information appended
#' @export
overlap_homology_and_result_dfs <- function(result_df, homology_df) {
  window_parts <- stringr::str_match(result_df$window, "^([^:]+):(\\d+)-(\\d+)$")

  # initialize granges object for result
  result_gr <- GenomicRanges::GRanges(
    seqnames = window_parts[, 2],
    ranges = IRanges::IRanges(
      start = as.integer(window_parts[, 3]),
      end = as.integer(window_parts[, 4])
    )
  )
  S4Vectors::mcols(result_gr) <- result_df

  # initialze granges object for homology df
  homology_gr <- GenomicRanges::GRanges(
    seqnames = homology_df$chromosome,
    ranges = IRanges::IRanges(
      start = homology_df$posit + 1L,
      width = homology_df$dna_width
    ),
    strand = homology_df$strand
  )
  S4Vectors::mcols(homology_gr) <- homology_df

  # find overlaps
  hits <- GenomicRanges::findOverlaps(
    query = result_gr,
    subject = homology_gr,
    ignore.strand = TRUE
  )

  # add columns to the result df related to homology
  result_df_new <- result_df |>
    dplyr::mutate(has_homology_hit = FALSE, homology_n_mismatches = NA_real_,
                  homology_n_bulges = NA_real_, homology_posit = NA_real_,
                  homology_strand = NA_character_, homology_dna = NA_character_,
                  homology_gRNA = NA_character_, homology_dna_width = NA_integer_)

  # create hit df, which records all the hits, alongside n_total_changes,
  if (length(hits) == 0L) {
    stop("No overlap between homology data frame and result data frame.")
  }

  hit_df <- data.frame(
    result_idx = S4Vectors::queryHits(hits),
    homology_idx = S4Vectors::subjectHits(hits)) |>
    dplyr::mutate(
      alignment_score = homology_df$n_mismatches[homology_idx] + 2 * homology_df$n_bulges[homology_idx]
    )

  n_hits <- tabulate(hit_df$result_idx, nbins = nrow(result_df_new))
  result_df_new$has_homology_hit <- n_hits > 0L

  best_hit_df <- hit_df |>
    dplyr::arrange(result_idx, alignment_score, homology_idx) |>
    dplyr::group_by(result_idx) |>
    dplyr::slice(1L) |>
    dplyr::ungroup()

  # combine information across data frames
  result_idx <- best_hit_df$result_idx
  homology_idx <- best_hit_df$homology_idx
  result_df_new$homology_n_mismatches[result_idx] <- homology_df$n_mismatches[homology_idx]
  result_df_new$homology_n_bulges[result_idx] <- homology_df$n_bulges[homology_idx]
  result_df_new$homology_posit[result_idx] <- homology_df$posit[homology_idx]
  result_df_new$homology_dna_width[result_idx] <- homology_df$dna_width[homology_idx]
  result_df_new$homology_strand[result_idx] <- homology_df$strand[homology_idx]
  result_df_new$homology_dna[result_idx] <- homology_df$dna[homology_idx]
  result_df_new$homology_gRNA[result_idx] <- homology_df$gRNA[homology_idx]

  # identify the cut site
  result_df_new_w_cut <- result_df_new |>
    dplyr::mutate(homology_cut_start = homology_posit + 1L + ifelse(homology_strand == "+", homology_dna_width - 4L, 3L),
                  homology_cut_end = homology_posit + 1L + ifelse(homology_strand == "+", homology_dna_width - 3L, 4L))

  return(result_df_new_w_cut)
}
