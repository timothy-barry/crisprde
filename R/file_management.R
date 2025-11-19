#' data_dir <- paste0(.get_config_path("LOCAL_CRISPR_DE_DATA_DIR"), "guideseq/count_tables/")
get_trt_cntrl_tab <- function(data_dir) {
  fs <- list.files(data_dir, full.names = FALSE)

  fs_mod <- gsub(pattern = "dsODN-only", replacement = "dsODNonly", x = fs) |>
    gsub(pattern = "ELANE-e1SD", replacement = "ELANEe1SD") |>
    gsub(pattern = "ELANE-e3SA", replacement = "ELANEe3SA") |>
    gsub(pattern = "-5uM", replacement = "")

  get_key_info_from_file <- function(file_name) {
    str_pieces <- strsplit(x = file_name, split = "-", fixed = TRUE)[[1]]
    cell_type <- str_pieces[1]
    crispr <- paste0(str_pieces[2], str_pieces[3])
    grna <- str_pieces[4]
    primer <- str_pieces[5]
    sample_id <- strsplit(x = str_pieces[6], split = "_", fixed = TRUE)[[1]][1]
    data.frame(cell_type = cell_type, crispr = crispr, grna = grna, primer = primer, sample_id = sample_id)
  }

  file_df <- lapply(fs_mod, FUN = get_key_info_from_file) |>
    data.table::rbindlist() |>
    dplyr::mutate(f = fs, f_mod = fs_mod)

  # get the treatment f; for each treatment f, determine its matched control f
  trt_sample_idxs <- which(file_df$grna != "dsODNonly")
  out <- matched_f_df <- lapply(trt_sample_idxs, FUN = function(trt_sample_idx) {
    curr_trt_info <- file_df[trt_sample_idx,c("cell_type", "crispr", "primer", "f")]
    cntrl_f <- dplyr::filter(file_df, cell_type == curr_trt_info$cell_type,
                             crispr == curr_trt_info$crispr,
                             primer == curr_trt_info$primer,
                             grna == "dsODNonly") |> dplyr::slice(1) |>
      dplyr::pull(f)
    data.frame(treat_f = curr_trt_info$f, cntrl_f = cntrl_f)
  }) |> data.table::rbindlist()

  return(out)
}
