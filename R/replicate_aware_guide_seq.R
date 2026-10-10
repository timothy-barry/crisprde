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


# helper function to convolve two pmfs; faster than the original implementation due to padding
convolve_pmfs_v2 <- function(a, b) {
  n <- length(a) + length(b) - 1L
  n_fft <- 2^ceiling(log2(n))

  a <- c(a, numeric(n_fft - length(a)))
  b <- c(b, numeric(n_fft - length(b)))

  out <- stats::fft(stats::fft(a) * stats::fft(b), inverse = TRUE)
  out <- Re(out[seq_len(n)]) / n_fft
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


fit_multirep_guideseq_occupancy <- function(Y_mat, incorporate_occupancy_info = TRUE) {
  MIN_NONZERO_COUNT <- 25L
  if (is.null(colnames(Y_mat))) warning("Y_mat must have column names (to identify the windows).")
  X <- Y_mat > 0
  storage.mode(X) <- "integer"

  nonzero_replicate_count <- rowSums(X)
  if (any(nonzero_replicate_count <= MIN_NONZERO_COUNT)) {
    offending_rows <- paste0(which(nonzero_replicate_count <= MIN_NONZERO_COUNT), collapse = ", ")
    msg <- paste0("Row ",  offending_rows, " has fewer than ", MIN_NONZERO_COUNT, " windows with a nonzero count. Consider dropping this sample or combining this sample with another (e.g., by pooling together primer chanels within a replicate).")
    warning(msg)
  }

  Omega <- as.matrix(generate_omega(nrow(Y_mat)))
  col_keys <- do.call(what = paste0, args = lapply(seq_len(nrow(X)), function(i) X[i, ]))
  names(col_keys) <- colnames(X)
  pi_hat <- NULL
  tbp_pattern_df <- NULL
  occupancy_pattern_map <- NULL

  if (incorporate_occupancy_info) {
    pi_hat <- fit_tbp_model(X)
    if (!is.null(pi_hat)) {
      tbp_pattern_df <- pmf_tbp(pi_hat)
    } else {
      warning("Cannot fit occupancy model; defaulting to count-only model.")
      incorporate_occupancy_info <- FALSE
    }
  }
  omega_keys <- apply(Omega, 1, paste0, collapse = "")
  occupancy_pattern_map <- match(col_keys, omega_keys)

  ret <- list(X = X,
              Omega = Omega,
              col_keys = col_keys,
              incorporate_occupancy_info = incorporate_occupancy_info,
              pi_hat = pi_hat,
              tbp_pattern_df = tbp_pattern_df,
              occupancy_pattern_map = occupancy_pattern_map)
  return(ret)
}


fit_multirep_guideseq_count_null <- function(Y_mat, c_grid) {
  # obtain the pilot fit for each row
  pilot_fit_list <- apply(X = Y_mat, MARGIN = 1, FUN = function(curr_row) {
    y_plus <- curr_row[curr_row > 0]
    y_plus_tab <- table(y_plus)
    y_compressed <- as.integer(names(y_plus_tab)) - 1L
    y_compressed_weights <- as.integer(y_plus_tab)
    shifted_fit_nb_pilot <- MASS::glm.nb(formula = y_compressed ~ 1, weights = y_compressed_weights)
    pilot <- c(mu = exp(shifted_fit_nb_pilot$coefficients[[1]]), theta = shifted_fit_nb_pilot$theta)
    l <- list(y_compressed = y_compressed, y_compressed_weights = y_compressed_weights, pilot = pilot)
  })
  # next, iterate over c_grid, producing the robust estimate for each c
  out <- lapply(X = c_grid, FUN = function(curr_c) {
    estimate_mat <- lapply(X = pilot_fit_list, FUN = function(curr_rep_pilot) {
      fit_fast <- fit_rob_nb_univariate(y = curr_rep_pilot$y_compressed,
                                        weights = curr_rep_pilot$y_compressed_weights,
                                        c.tukey.beta = curr_c,
                                        c.tukey.sigma = curr_c,
                                        pilot = curr_rep_pilot$pilot)
    }) |> dplyr::bind_rows() |> as.matrix()
    rownames(estimate_mat) <- rownames(Y_mat)
    return(estimate_mat)
  }) |> setNames(c_grid)
  return(out)
}


get_p_values_given_test_stats_prob_vector <- function(test_stat_v_in, right_tail_prob_v_in) {
  p_vals <- numeric(length(test_stat_v_in))
  neg_idx <- which(test_stat_v_in < 0)
  ok_idx <- which(test_stat_v_in >= 0)
  p_vals[neg_idx] <- 1
  tail_idx <- as.integer(test_stat_v_in[ok_idx]) + 1L
  too_large <- tail_idx > length(right_tail_prob_v_in)
  p_vals[ok_idx[too_large]] <- right_tail_prob_v_in[length(right_tail_prob_v_in)]
  p_vals[ok_idx[!too_large]] <- right_tail_prob_v_in[tail_idx[!too_large]]
  return(p_vals)
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
#' - `window` (a string indicating the cluster to which a given base belongs)
#' - `cluster_chr` (chromosome of the group)
#' - `min_cluster_coord` (minimum coordinate of the group)
#' - `max_cluster_coord` (maximum coordinate of the group)
#' @export
#'
#' @examples
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' count_df <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#'  dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & treated & replicate_id %in% 1:2) |>
#'  dplyr::filter(chr != "chrM") |>
#'  dplyr::select(chr, coord, strand, umi_count, primer_type, replicate_id)
#' clustered_count_df <- cluster_loci(count_df)
cluster_loci <- function(count_df, thresh = 100L, padding = 5L) {
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
    dplyr::mutate(cluster_chr = chr[1],
                  min_cluster_coord = min(coord) - padding,
                  max_cluster_coord = max(coord) + padding) |>
    dplyr::ungroup() |>
    dplyr::mutate(window = paste0(cluster_chr, ":", min_cluster_coord, "-", max_cluster_coord),
                  group_id = NULL, simple_dist = NULL)
  return(count_df_w_dist_and_group_string)
}


#' Construct replicate count table
#'
#' Takes the output of `cluster_loci()` as input; outputs a count matrix for statistical modeling
#'
#' @param clustered_count_df output of `cluster_loci()`, with columns `window`, `umi_count`, and the columns specified by `channel_axes`.
#' @param channel_axes columns whose unique tuples define evidence channels.
#'
#' @returns an integer matrix with evidence channels in the rows and windows in the columns. An entry of the matrix indicates the number of UMIs observed within a given window and channel.
#' @export
#'
#' @examples
#' elane_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_elane/")
#' count_df <- readRDS(paste0(elane_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#'  dplyr::filter(cell_type == "CD34" & cas9_variant == "wt_cas9" & treated & replicate_id %in% 1:2, chr != "chrM") |>
#'  dplyr::select(chr, coord, strand, umi_count, primer_type, replicate_id)
#' clustered_count_df <- cluster_loci(count_df)
#' Y_mat <- construct_replicate_count_table(clustered_count_df)
construct_replicate_count_table <- function(clustered_count_df, channel_axes = c("replicate_id", "primer_type")) {
  x <- lapply(X = channel_axes, FUN = function(col_name) {
    clustered_count_df[[col_name]] |> as.character()
  })
  clustered_count_df$channel_id <- Reduce(f = function(a, b) paste0(a, "_", b), x = x)

  # sum over UMIs within a given (window, channel) pair
  collapsed_count_df <- clustered_count_df |>
    dplyr::group_by(channel_id, window) |>
    dplyr::summarize(umi_count = sum(umi_count), .groups = "drop") |>
    dplyr::arrange(channel_id, window)

  # construct Y_mat
  replicate_idx <- match(x = collapsed_count_df$channel_id, unique(collapsed_count_df$channel_id))
  group_idx <- match(x = collapsed_count_df$window, unique(collapsed_count_df$window))
  Y_mat <- Matrix::sparseMatrix(i = replicate_idx,
                                j = group_idx,
                                x = collapsed_count_df$umi_count) |>
    as.matrix()
  rownames(Y_mat) <- unique(collapsed_count_df$channel_id)
  colnames(Y_mat) <- unique(collapsed_count_df$window)
  return(Y_mat)
}


#' Load CRISPRitz output
#'
#' @param targets_file_path file path to the targets.txt file outputted by CRISPRitz
#'
#' @returns a data frame containing the CRISPRitz target output
#' @export
#'
#' @examples
#' homology_df <- load_crispritz_output("/Users/timbarry/research_offsite/external/bauer-lab/guideseq_elane/crispritz_CCCCGGCAGAAACGTCCGCG.hg38.targets.txt")
load_crispritz_output <- function(targets_file_path) {
  df <- readr::read_delim(file = targets_file_path) |>
    dplyr::rename("bulge_type" = "#Bulge type", "n_mismatches" = "Mismatches",
                  "n_bulges" = "Bulge Size", "n_total_changes" = "Total",
                  "gRNA" = "crRNA", "dna" = "DNA", "chromosome" = "Chromosome",
                  "posit" = "Position", "cluster_posit" = "Cluster Position",
                  "strand" = "Direction")
  protospacer_width <- nchar(gsub(pattern = "-", replacement = "", x = df$dna)) - 3L
  df$protospacer_width <- protospacer_width
  return(df)
}


#' Load N run bed file
#'
#' Load bed file storing the positions of N-runs in the reference genome.
#'
#' @param n_run_bed_file_path path to N-run bed file
#'
#' @returns a data frame storing the positions of the N-runs
#' @export
#'
#' @examples
#' n_run_bed_file_path <- "/Users/timbarry/research_offsite/ref_genome_dir/hg38_N_runs_min10.bed"
#' n_run_df <- load_n_run_bed(n_run_bed_file_path)
load_n_run_bed <- function(n_run_bed_file_path) {
  df <- readr::read_delim(file = n_run_bed_file_path, col_names = FALSE)
  colnames(df) <- c("chromosome", "start", "end", "feature", "score", "strand")
  df |> dplyr::mutate(start = start + 1L)
}


#' Load ENCODE blacklist bed
#'
#' @param encode_blacklist_file_path filepath to ENCODE blacklist
#'
#' @returns a data frame storing the positions of the blacklisted regions
#' @export
#'
#' @examples
#' encode_blacklist_file_path <- "/Users/timbarry/research_offsite/ref_genome_dir/hg38-blacklist.v2.bed"
load_encode_blacklist_bed <- function(encode_blacklist_file_path) {
  df <- readr::read_delim(file = encode_blacklist_file_path, col_names = FALSE)
  colnames(df) <- c("chromosome", "start", "end", "classification")
  df |> dplyr::mutate(start = start + 1L)
}

#' Annotate clustered count df with homology
#'
#' Add columns:
#' - `homology_has_hit` (indicating whether there is a homology hit inside the window)
#' - `overlaps_n_run` (indicating whether the window overlaps an N-run in the reference genome)
#' - `overlaps_encode_blacklist` (indicating whether the window overlaps an ENCODE blacklist region)
#' - `homology_n_mismatches` (indicating the number of mismatches between aligned spacer and protospacer sequence)
#' - `homology_n_bulges` (indicating number of bulges between aligned spacer and protospacer sequence)
#' - `homology_bulge_type` (indicating the CRISPRitz bulge type)
#' - `homology_cfd` (indicating the CFD-like homology score)
#' - `homology_chromosome` (chromosome of the aligned protospacer)
#' - `homology_strand` (indicating whether the protospacer is on the plus or minus strand)
#' - `homology_dna` (aligned protospacer sequence)
#' - `homology_gRNA` (aligned spacer sequence)
#' - `homology_cut_start` and `homology_cut_end` (predicted cut-site coordinates)
#' - `homology_modal_base_cut_distance` (distance between the modal base and predicted cut site)
#' - `homology_alignment_score` (combined homology and cut-site-distance score)
#'
#' Windows without a CRISPRitz hit have CFD zero, distance Inf, and alignment score zero.
#'
#' @param clustered_count_df output of `cluster_loci()`
#' @param homology_df optional output of `load_crispritz_output()`; if supplied, windows are annotated for overlap with CRISPRitz hits
#' @param n_run_df optional output of `load_n_run_bed()`; if supplied, windows are annotated for overlap with N-runs
#' @param encode_blacklist_df optional output of `load_encode_blacklist_bed()`; if supplied, windows are annotated for overlap with ENCODE blacklist regions
#' @param gamma exponential distance-decay coefficient for ranking candidate alignments
#'
#' @examples
#' homology_df <- load_crispritz_output("/Users/timbarry/research_offsite/external/bauer-lab/guideseq_bcl11a/1620_crispritz_spRY_bcl11a_windows.hg38.targets.txt")
#' n_run_df <- load_n_run_bed("/Users/timbarry/research_offsite/ref_genome_dir/hg38_N_runs_min10.bed")
#' encode_blacklist_df <- load_encode_blacklist_bed("/Users/timbarry/research_offsite/ref_genome_dir/hg38-blacklist.v2.bed")
#' bcl11a_dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_bcl11a/")
#' clustered_count_df <- readRDS(paste0(bcl11a_dir, "count_tables_no_multimap/combined_count_df.rds")) |>
#' dplyr::filter(grna == "1620" & treated & chr != "chrM") |>
#' dplyr::select(chr, coord, strand, umi_count, primer_type, replicate_id) |>
#' cluster_loci()
#' annotated_clustered_count_df <- annotate_clustered_count_df(clustered_count_df = clustered_count_df,
#'   homology_df = homology_df, n_run_df = n_run_df, encode_blacklist_df = encode_blacklist_df)
#' @export
annotate_clustered_count_df <- function(clustered_count_df, homology_df = NULL, n_run_df = NULL,
                                        encode_blacklist_df = NULL, gamma = log(20)/7) {
  # compute window df by computing a summary over clustered_count_df
  window_df <- clustered_count_df |>
    dplyr::group_by(window, coord) |>
    dplyr::summarize(umi_count_per_base = sum(umi_count),
                     chr = chr[1],
                     max_cluster_coord = max_cluster_coord[1],
                     min_cluster_coord = min_cluster_coord[1]) |>
    dplyr::summarize(modal_base = coord[umi_count_per_base == max(umi_count_per_base)][1],
                     window_n_occupied_bases = length(unique(coord)),
                     window_width = max_cluster_coord[1] - min_cluster_coord[1] + 1L,
                     chr = chr[1],
                     min_cluster_coord = min_cluster_coord[1],
                     max_cluster_coord = max_cluster_coord[1])

  # construct window gr
  window_gr <- GenomicRanges::GRanges(
    seqnames = window_df$chr,
    ranges = IRanges::IRanges(start = window_df$min_cluster_coord, end = window_df$max_cluster_coord)
  )

  append_simple_overlap <- function(window_df, window_gr, df_in, colname) {
    gr_in <- GenomicRanges::GRanges(
      seqnames = df_in$chromosome,
      ranges = IRanges::IRanges(start = df_in$start, end = df_in$end)
    )
    seqlevels <- union(GenomeInfoDb::seqlevels(window_gr), GenomeInfoDb::seqlevels(gr_in))
    GenomeInfoDb::seqlevels(window_gr) <- seqlevels
    GenomeInfoDb::seqlevels(gr_in) <- seqlevels
    hits <- GenomicRanges::findOverlaps(query = window_gr, subject = gr_in,
                                        ignore.strand = TRUE)
    window_df[[colname]] <- FALSE
    window_df[[colname]][S4Vectors::queryHits(hits)] <- TRUE
    return(window_df)
  }
  if (!is.null(n_run_df)) {
    window_df <- append_simple_overlap(window_df = window_df,
                                       window_gr = window_gr,
                                       df_in = n_run_df,
                                       colname = "overlaps_n_run")
  }
  if (!is.null(encode_blacklist_df)) {
    window_df <- append_simple_overlap(window_df = window_df,
                                       window_gr = window_gr,
                                       df_in = encode_blacklist_df,
                                       colname = "overlaps_encode_blacklist")
  }

  # add homology information if homology_df provided
  if (!is.null(homology_df)) {
    # construct homology gr object
    homology_gr <- GenomicRanges::GRanges(
      seqnames = homology_df$chromosome,
      ranges = IRanges::IRanges(
        start = homology_df$posit + 1L + ifelse(homology_df$strand == "-", 3L, 0L),
        width = homology_df$protospacer_width),
      strand = homology_df$strand)

    # find overlaps between homology gr and window gr; retain alignments overlapping a GUIDE-seq window; append cfd score
    crispritz_hits <- GenomicRanges::findOverlaps(query = window_gr,
                                                  subject = homology_gr,
                                                  ignore.strand = TRUE)
    if (length(crispritz_hits) == 0L) {
      window_df <- window_df |>
        dplyr::mutate(homology_has_hit = FALSE,
                      homology_bulge_type = NA_character_, homology_gRNA = NA_character_,
                      homology_dna = NA_character_, homology_chromosome = NA_character_,
                      homology_strand = NA_character_, homology_n_mismatches = NA_real_,
                      homology_n_bulges = NA_real_, homology_cut_start = NA_real_,
                      homology_cut_end = NA_real_, homology_cfd = NA_real_,
                      homology_modal_base_cut_distance = NA_real_, homology_alignment_score = NA_real_)
    } else {
      crispritz_subject_hits <- S4Vectors::subjectHits(crispritz_hits)
      crispritz_query_hits <- S4Vectors::queryHits(crispritz_hits)
      homology_df_sub <- homology_df[crispritz_subject_hits,] |>
        dplyr::mutate(homology_cut_start = posit + 1L + ifelse(strand == "+", protospacer_width - 4L, 5L),
                      homology_cut_end = posit + 1L + ifelse(strand == "+", protospacer_width - 3L, 6L)) |>
        dplyr::mutate(homology_cfd = calculate_cfd_score(homology_dna = dna, homology_gRNA = gRNA, homology_has_hit = TRUE))

      # for each window (containing a CRISPRitz hit), score each alignment and return the best
      window_idxs_with_hit <- unique(crispritz_query_hits)
      homology_df <- lapply(X = window_idxs_with_hit, FUN = function(i) {
        curr_window <- window_df[i,]
        curr_crispritz_candidates <- homology_df_sub[crispritz_query_hits == i,] |>
          dplyr::mutate(homology_modal_base_cut_distance = pmin(abs(curr_window$modal_base - homology_cut_start),
                                                                abs(curr_window$modal_base - homology_cut_end))) |>
          dplyr::mutate(homology_alignment_score = compute_alignment_scores(cfds = homology_cfd, distances = homology_modal_base_cut_distance, gamma = gamma))
        # find the alignment with the best score, tie-breaking by distance
        best_alignment <- curr_crispritz_candidates |>
          dplyr::arrange(dplyr::desc(homology_alignment_score), homology_modal_base_cut_distance, dplyr::desc(homology_cfd)) |>
          dplyr::slice(1L) |>
          dplyr::select(bulge_type, gRNA, dna, chromosome, strand, n_mismatches, n_bulges,
                        homology_cut_start, homology_cut_end, homology_cfd, homology_modal_base_cut_distance,
                        homology_alignment_score) |>
          dplyr::rename(homology_bulge_type = bulge_type, homology_gRNA = gRNA, homology_dna = dna,
                        homology_chromosome = chromosome, homology_strand = strand, homology_n_mismatches = n_mismatches,
                        homology_n_bulges = n_bulges) |>
          dplyr::mutate(window = curr_window$window)
      }) |> data.table::rbindlist()
      window_df <- dplyr::left_join(window_df, homology_df, by = "window")
    }
  }

  # prepare final output
  window_df <- window_df |> dplyr::select(-min_cluster_coord, -max_cluster_coord, -chr)
  clustered_count_df <- dplyr::left_join(x = clustered_count_df, y = window_df, by = "window")
  if (!is.null(homology_df)) {
    clustered_count_df <- clustered_count_df |>
      dplyr::mutate(homology_has_hit = ifelse(is.na(homology_bulge_type), FALSE, TRUE),
                    homology_cfd = ifelse(homology_has_hit, homology_cfd, 0),
                    homology_modal_base_cut_distance = ifelse(homology_has_hit, homology_modal_base_cut_distance, Inf),
                    homology_alignment_score = ifelse(homology_has_hit, homology_alignment_score, 0))
  }

  return(clustered_count_df)
}


remove_duplicate_umis_with_shared_base <- function(count_df) {
  count_df$row_id <- seq_len(nrow(count_df))
  duplicated_umi_df <- data.frame(row_id = rep(count_df$row_id, lengths(count_df$umis)),
                                  chr = rep(count_df$chr, lengths(count_df$umis)),
                                  coord = rep(count_df$coord, lengths(count_df$umis)),
                                  umi = unlist(count_df$umis, use.names = FALSE)) |>
    dplyr::group_by(chr, coord) |>
    dplyr::mutate(duplicated = duplicated(umi) | duplicated(umi, fromLast = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::filter(duplicated)

  # a list of duplicated umis, ordered by row
  duplicated_umis_by_row <- split(x = duplicated_umi_df$umi, f = duplicated_umi_df$row_id)

  # remove the duplicated umis
  for (i in seq_along(duplicated_umis_by_row)) {
    umi_to_remove <- duplicated_umis_by_row[[i]]
    row_idx <- names(duplicated_umis_by_row[i]) |> as.integer()
    curr_umi_vector <- count_df[[row_idx, "umis"]]
    keep_v <- !(curr_umi_vector %in% umi_to_remove)

    # update n_umis, umis, and reads per umi
    count_df[[row_idx, "umi_count"]] <- sum(keep_v)
    if (count_df[[row_idx, "umi_count"]] >= 1L) {
      filtered_umi_vector <- curr_umi_vector[keep_v]
      filtered_n_reads_per_umi <- count_df[[row_idx, "n_reads_per_umi"]][keep_v]
      count_df[row_idx, umis := I(filtered_umi_vector)]
      count_df[row_idx, n_reads_per_umi := I(filtered_n_reads_per_umi)]
      count_df[row_idx, total_read_count := sum(filtered_n_reads_per_umi)]
    }
  }
  count_df <- count_df |> dplyr::filter(umi_count >= 1L)
  count_df$row_id <- NULL
  return(count_df)
}
