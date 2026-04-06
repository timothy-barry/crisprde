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
#' data_list <- readRDS("/Users/timbarry/research_offsite/external/bauer-lab/rhampseq_bcl11a/multivariate_count_tables/combined/CRISPResso_on_1617_OT_0068.rds")
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
#' data_list <- readRDS("/Users/timbarry/research_offsite/external/bauer-lab/rhampseq_bcl11a/multivariate_count_tables/combined/CRISPResso_on_1617_OT_0068.rds")
#' count_matrix <- data_list$count_matrix
#' allele_df <- data_list$allele_df
#' covariate_df <- data_list$covariate_df
#' collapsed_matrix <- collapse_count_matrix_by_allele(count_matrix, allele_df)
#' collapsed_res <- compute_basic_statistics_on_multivariate_amplicon_seq_data(collapsed_matrix, covariate_df)
#'
collapse_count_matrix_by_allele <- function(count_matrix, allele_df,
                                            bucket_breaks = c(0L, 1L, 2L, 3L, 4L, 5L, 7L, 9L, Inf),
                                            bucket_labels = c("one", "two", "three", "four", "five", "six-seven", "eight-nine", "10+")) {
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
  ret <- cbind(unmutated_count_vect, count_matrix)
  return(ret)
}
