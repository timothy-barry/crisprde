parse_crispresso_vector <- function(position_strings) {
  lapply(X = strsplit(gsub("\\[|\\]", "", position_strings), ","),
         FUN = function(position_values) as.integer(position_values) + 1L)
}

has_substitution_helper <- function(allele_table, substitution_window_reference_positions, reference_amplicon, from_nuc, to_nuc) {
  # parse substitution positions and alignment-column mappings into 1-based coordinates
  substitution_reference_positions_by_allele <- parse_crispresso_vector(allele_table$all_substitution_positions)
  alignment_to_reference_positions_by_allele <- parse_crispresso_vector(allele_table$ref_positions)
  aligned_read_bases_by_allele <- allele_table$Aligned_Sequence |> strsplit(split = "")

  # flag alleles with at least one specified substitution in the window
  has_specified_substitution <- vapply(X = seq_along(substitution_reference_positions_by_allele), FUN = function(allele_index) {
    substitution_reference_positions <- substitution_reference_positions_by_allele[[allele_index]]
    if (length(substitution_reference_positions) == 0L) {
      FALSE
    } else {
      # prepare to index into the aligned sequence in the substitution window
      alignment_to_reference_positions <- alignment_to_reference_positions_by_allele[[allele_index]]
      substitution_alignment_columns <- match(x = substitution_reference_positions, table = alignment_to_reference_positions)
      aligned_read_bases <- aligned_read_bases_by_allele[[allele_index]]
      any(substitution_reference_positions %in% substitution_window_reference_positions &
            reference_amplicon[substitution_reference_positions] == from_nuc &
            aligned_read_bases[substitution_alignment_columns] == to_nuc)
    }
  }, FUN.VALUE = logical(1L))
}

reverse_complement <- function(bases) chartr("ACGT", "TGCA", rev(bases))

#' Compute base editing count table
#'
#' This function takes a CRISPResso directory as input. It outputs a data frame with one row containing the following columns:
#'  - n_reads: total number of (aligned) reads for this sample.
#'  - n_reads_w_intended_substitution: number of reads harboring the intended substitution.
#'  - n_reads_w_bystander_substitution: number of reads harboring a bystander substitution.
#'  - n_reads_w_indel: number of reads harboring an indel.
#'  - n_reads_w_modification: number of reads harboring an intended substitution OR bystander substitution OR indel.
#'
#' Identify modification locations using 1-based indexing relative to the distal end of the protospacer sequence. For example, the default `indel_positions` is `c(17L, 18L)` (i.e., the bases immediately adjacent to the nicking site), while the default `bystander_substitution_positions` is `seq(3L, 11)` (i.e., positions 3-11 from the distal end of the protospacer).
#'
#' Specify the substitution type using the `from_nuc` and `to_nuc` arguments; for an adenine base editor, `from_nuc = "A"` and `to_nuc = "G"`; for a cytosine base editor, `from_nuc = "C"` and `to_nuc = "T"`.
#'
#'  Note that `n_reads_w_modification` is not simply the sum over `n_reads_w_intended_substitution`, `n_reads_w_bystander_substitution`, and `n_reads_w_indel`, as some reads may harbor multiple modifications of different types (e.g., an intended substitution and an indel).
#'
#' @param crispresso_dir file path to a CRISPResso output directory
#' @param indel_positions bases over which to quantify indels; default bases flanking PAM
#' @param bystander_substitution_positions bases over which to quantify bystander substitutions; default positions 3-11 from distal end of protospacer
#' @param intended_substitution_positions bases over which to quantify intended substitutions; default integer(0L) (i.e., empty set)
#' @param from_nuc original nucleotide; default "A"
#' @param to_nuc converted nucleotide; default "G"
#'
#' @returns a 1-row data frame with columns `n_reads`, `n_reads_w_intended_substitution`, `n_reads_w_bystander_substitution`, `n_reads_w_indel`, `n_reads_w_modification`
#' @export
#'
#' @examples
#' crispresso_dir <- "/Users/timbarry/research_offsite/external/scratch/novaseq_12_elane_mrna_dose_response/chr19_853252_853275/crispresso_outputs/CRISPResso_on_Jing_i507_i733_on_target_A8toG/"
#' count_df <- compute_base_editing_count_table(crispresso_dir = crispresso_dir,
#'                                              bystander_substitution_positions = 10L,
#'                                              intended_substitution_positions = 8L)
compute_base_editing_count_table <- function(crispresso_dir, indel_positions = c(17L, 18L),
                                             bystander_substitution_positions = seq(3L, 11),
                                             intended_substitution_positions = integer(0L),
                                             from_nuc = "A", to_nuc = "G") {
  # extract reference amplicon, protospacer, and aligned-read metadata
  crispresso_metadata <- jsonlite::read_json(file.path(crispresso_dir, "CRISPResso2_info.json"))
  amplicon_metadata <- crispresso_metadata$results$refs[[1]]
  reference_amplicon <- (amplicon_metadata$sequence |> strsplit(split = ""))[[1]]
  protospacer_seq <- (amplicon_metadata$sgRNA_sequences[[1]] |> strsplit(split = ""))[[1]]
  protospacer_interval_0based <- unlist(amplicon_metadata$sgRNA_intervals[[1]])
  protospacer_interval_1based <- protospacer_interval_0based + 1L
  n_aligned_reads <- crispresso_metadata$results$alignment_stats$counts_total[[1]]

  # determine whether protospacer is on plus or minus strand
  protospacer_on_amplicon <- reference_amplicon[protospacer_interval_1based[1]:protospacer_interval_1based[2]]
  if (all(protospacer_seq == protospacer_on_amplicon)) {
    start <- protospacer_interval_0based[2] - 19L
    intended_substitution_window_reference_positions <- start + intended_substitution_positions
    bystander_substitution_window_reference_positions <- start + bystander_substitution_positions
    indel_window_reference_positions <- start + indel_positions
  } else if (all(protospacer_seq == reverse_complement(protospacer_on_amplicon))) {
    start <- protospacer_interval_1based[1]
    intended_substitution_window_reference_positions <- start + 20L - intended_substitution_positions
    bystander_substitution_window_reference_positions <- start + 20L - bystander_substitution_positions
    indel_window_reference_positions <- start + 20L - indel_positions
    from_nuc <- reverse_complement(from_nuc)
    to_nuc <- reverse_complement(to_nuc)
  } else {
    stop("Protospacer sequence not recognized on amplicon.")
  }

  # load the allele frequency table
  allele_table <- paste0(crispresso_dir, "Alleles_frequency_table.zip") |>
    data.table::fread() |> dplyr::rename("n_reads" = "#Reads")
  has_bystander_substitution <- has_substitution_helper(allele_table = allele_table,
                                                        substitution_window_reference_positions = bystander_substitution_window_reference_positions,
                                                        reference_amplicon = reference_amplicon, from_nuc = from_nuc, to_nuc = to_nuc)
  has_intended_substitution <- has_substitution_helper(allele_table = allele_table,
                                                       substitution_window_reference_positions = intended_substitution_window_reference_positions,
                                                       reference_amplicon = reference_amplicon, from_nuc = from_nuc, to_nuc = to_nuc)

  # flag alleles with an indel overlapping the cut site
  deletion_reference_positions_by_allele <- parse_crispresso_vector(allele_table$all_deletion_positions)
  has_deletion <- vapply(X = deletion_reference_positions_by_allele, FUN = function(deletion_position) {
    any(deletion_position %in% indel_window_reference_positions)
  }, FUN.VALUE = logical(1L))

  # flag alleles with insertion overlapping the cut site
  left_insertion_reference_positions_by_allele <- parse_crispresso_vector(allele_table$all_insertion_left_positions)
  has_insertion <- vapply(X = left_insertion_reference_positions_by_allele, FUN = function(insertion_left_position) {
    any(insertion_left_position %in% indel_window_reference_positions &
          (insertion_left_position + 1L) %in% indel_window_reference_positions
    )
  }, FUN.VALUE = logical(1L))
  has_indel <- has_deletion | has_insertion

  # finally, compute number of reads with a specified substitution, indel, or either
  ret <- allele_table |>
    dplyr::mutate(has_intended_substitution = has_intended_substitution,
                  has_bystander_substitution = has_bystander_substitution,
                  has_indel = has_indel,
                  has_modification = has_intended_substitution | has_bystander_substitution | has_indel) |>
    dplyr::summarize(n_reads_w_intended_substitution = sum(n_reads[has_intended_substitution]),
                     n_reads_w_bystander_substitution = sum(n_reads[has_bystander_substitution]),
                     n_reads_w_indel = sum(n_reads[has_indel]),
                     n_reads_w_modification = sum(n_reads[has_modification])) |>
    dplyr::mutate(n_reads = n_aligned_reads) |> dplyr::relocate(n_reads)
  return(ret)
}
