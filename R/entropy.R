compute_site_entropy <- function(count_df, max_chain_link = 25, chrs_to_keep = seq(1L, 22L)) {
  count_df <- count_df |> dplyr::filter(chr %in% chrs_to_keep)
  clustered_count_df <- compute_distances_between_occupied_loci(count_df = count_df, thresh = max_chain_link)
  entropy_df <- clustered_count_df |>
    dplyr::group_by(group_id) |>
    dplyr::do(compute_entropy(.))
}

compute_entropy <- function(df) {
  # df <- clustered_count_df |> dplyr::filter(group_id == 28)
  count <- df$count
  p <- count/sum(count)
  K <- length(p)
  data.frame(entropy = -sum(p * log(p))/K)
}
