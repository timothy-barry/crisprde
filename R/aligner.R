#' Align a spacer sequence to the genome
#'
#' Aligns a short (20 bp) spacer sequence to the genome.
#'
#' @param grna_spacer gRNA spacer sequence (written in 5' -> 3' direction)
#' @param pam (optional) the PAM sequence
#' @param ref_genome file path to the reference genome
#'
#' @returns a GAlignments object containing the alignment
#' @export
#'
#' @examples
#' alignment <- align_spacer_to_genome("ACTGATAGGGGTCGCGGTAG")
align_spacer_to_genome <- function(grna_spacer, ref_genome = paste0(.get_config_path("REF_GENOME_DIR"), "hg38_main_chroms"), pam = NA) {
  if (!is.na(pam) && "N" %in% strsplit(x = pam, split = "")[[1]]) pam <- NA
  if (methods::is(grna_spacer, "DNAString")) grna_spacer <- as.character(grna_spacer)
  if (!is.na(pam)) grna_spacer <- paste0(grna_spacer, pam)

  temp_dir <- tempdir()
  temp_fasta <- paste0(temp_dir, "/spacer.fa")
  temp_sam <- paste0(temp_dir, "/spacer.sam")
  writeLines(c("@spacer", grna_spacer, "+", paste0(rep("I", nchar(grna_spacer)), collapse = "")), temp_fasta)
  bowtie_command <- paste0("cd ", temp_dir, "; bowtie2 -x ", ref_genome, " -U spacer.fa --score-min L,0,0 --mp 1000,1000 -k 100 -S spacer.sam")
  system(bowtie_command)
  bowtie_command_2 <- paste0("cd ", temp_dir, "; samtools view -bS spacer.sam > spacer.bam; samtools sort spacer.bam -o spacer.sorted.bam; samtools index spacer.sorted.bam")
  system(bowtie_command_2)

  # read SAM file
  alignment <- GenomicAlignments::readGAlignments(file = paste0(temp_dir, "/spacer.sorted.bam"))
  return(alignment)
}


#' Write a mini reference genome from count-table windows
#'
#' Clusters occupied bases, pads each window, and writes one FASTA file per
#' window. Padded coordinates are clipped to chromosome boundaries and encoded
#' as `chr:start-end` in filenames and sequence headers (1-based, inclusive).
#' Windows are kept separate even when their padded ranges overlap.
#'
#' @param count_df Count table containing `chr` and `coord` columns.
#' @param reference_fasta Path to an indexed reference FASTA file. Defaults to
#'   `hg38_main_chroms.fa` in `.get_config_path("REF_GENOME_DIR")`.
#' @param output_dir Directory for the mini reference genome. Existing `.fa`
#'   files in this directory are replaced.
#' @param thresh Maximum distance between occupied bases in the same cluster.
#' @param padding Number of bases to add on each side of a window.
#' @param exclude_chr_M Whether to exclude rows whose chromosome is `chrM`.
#'
#' @returns Invisibly returns `NULL`.
#' @export
write_mini_reference_genome <- function(count_df, reference_fasta = file.path(.get_config_path("REF_GENOME_DIR"), "hg38_main_chroms.fa"),
                                       output_dir, thresh = 100L, padding = 100L, exclude_chr_M = TRUE) {
  if (exclude_chr_M) count_df <- dplyr::filter(count_df, chr != "chrM")
  window_df <- cluster_loci(count_df, thresh = thresh, padding = 0L) |>
    dplyr::distinct(cluster_chr, min_cluster_coord, max_cluster_coord)

  fa <- Rsamtools::FaFile(reference_fasta)
  reference_index <- Rsamtools::scanFaIndex(fa)
  chromosome_lengths <- stats::setNames(GenomicRanges::width(reference_index),
                                        as.character(GenomicRanges::seqnames(reference_index)))
  windows <- GenomicRanges::GRanges(
    seqnames = window_df$cluster_chr,
    ranges = IRanges::IRanges(
      start = pmax(1L, window_df$min_cluster_coord - padding),
      end = pmin(chromosome_lengths[window_df$cluster_chr],
                 window_df$max_cluster_coord + padding)))
  window_seq <- Biostrings::getSeq(fa, windows)
  names(window_seq) <- paste0(GenomicRanges::seqnames(windows), ":",
                             GenomicRanges::start(windows), "-", GenomicRanges::end(windows))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  file.remove(list.files(output_dir, pattern = "\\.fa$", full.names = TRUE))
  for (i in seq_along(window_seq)) {
    Biostrings::writeXStringSet(window_seq[i],
                               filepath = file.path(output_dir, paste0(names(window_seq)[i], ".fa")))
  }
  invisible(NULL)
}


#' Map CRISPRitz output back to reference-genome coordinates
#'
#' Reads CRISPRitz output for a mini reference genome written by
#' `write_mini_reference_genome()`. Extracts the chromosome and 1-based window
#' start from each `chr:start-end` sequence name and adds `start - 1` to both
#' `Position` and `Cluster Position`, preserving CRISPRitz's 0-based convention.
#' All other columns and rows are retained without guide-specific processing.
#'
#' @param targets_file_path Path to the CRISPRitz targets TSV file.
#' @param output_file_path Path at which to write the remapped targets TSV file.
#'
#' @returns Invisibly returns `NULL`.
#' @export
remap_crispritz_output <- function(targets_file_path, output_file_path) {
  df <- readr::read_tsv(targets_file_path, show_col_types = FALSE)
  window_start <- as.integer(sub(".*:([0-9]+)-[0-9]+$", "\\1", df$Chromosome))
  df$Chromosome <- sub(":[0-9]+-[0-9]+$", "", df$Chromosome)
  df$Position <- df$Position + window_start - 1L
  df$`Cluster Position` <- df$`Cluster Position` + window_start - 1L
  readr::write_tsv(df, output_file_path)
  invisible(NULL)
}
