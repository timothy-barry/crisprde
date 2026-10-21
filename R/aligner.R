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
