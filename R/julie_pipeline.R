#' @importFrom stats setNames
#' @importFrom methods as
#' @rawNamespace import(S4Vectors, except=c(fold, values, rename))
#' @rawNamespace import(IRanges, except=values)
#' @rawNamespace import(BSgenome, except=export)
#' @rawNamespace import(GenomicRanges, except=values)
#' @rawNamespace import(BiocGenerics, except=c(var, sd))
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics hist par
#' @importFrom GenomicAlignments readGAlignmentsList first
#' last readGAlignmentPairs
#' @importFrom data.table fread
#' @importFrom parallel makeCluster stopCluster detectCores
#' parLapply
#' @importFrom Rsamtools ScanBamParam BamFile bamFlagTest
#' @importFrom tools file_ext
#' @importFrom dplyr select mutate add_count filter slice_sample group_by count '%>%'
#' @examples
#' alignment.inputfile <- "/Users/timbarry/research_offsite/external/bauer-lab/guideseq/bam_files/293T-SpRY-Cas9-1620-GSPneg-S79_S87_L001_alignment.bam"
#' umi.inputfile <- "/Users/timbarry/research_offsite/external/bauer-lab/guideseq/umi_tables/293T-SpRY-Cas9-1620-GSPneg-S79_S87_L001_umi_table.txt"
getUniqueCleavageEvents <-
  function(alignment.inputfile,
           umi.inputfile,
           max.paired.distance = 1000,
           max.R1.len = 130,
           max.R2.len = 130,
           same.chromosome = TRUE,
           distance.inter.chrom = -1,
           min.R1.mapped = 20,
           min.R2.mapped = 20,
           apply.both.min.mapped = FALSE,
           max.duplicate.distance = 0L,
           umi.plus.R1start.unique = TRUE,
           umi.plus.R2start.unique = TRUE,
           min.umi.count = 5L,
           max.umi.count = 100000L,
           min.read.coverage = 1L,
           n.cores.max = 6,
           removeDuplicate = TRUE,
           ignoreTagmSite = FALSE,
           ignoreUMI = FALSE) {
    library(GenomicRanges)
    library(BSgenome)
    library(BiocGenerics)
    library(S4Vectors)
    library(IRanges)

    # load UMI table
    umi <- as.data.frame(data.table::fread(umi.inputfile, sep = "\t", colClasses = "character", header = FALSE))
    # align BAM alignments via importBAMAlignments (below)
    align <- importBAMAlignments(alignment.inputfile)
    # set umi column names
    colnames(umi) <- c("readName", "UMI")
    # remote @ symbol from readname
    umi$readName <- gsub("^@", "", umi$readName)
    # keep only reads whose width is less than 130 (if some reads do not contain the dsODN tag, they are filtered here)
    align <- subset(align, (width.first > 0L & width.first <= max.R1.len) | (width.last > 0L & width.last <= max.R2.len))
    # add umi to align
    align.umi <- merge(align, umi)
    # remove reads with UMIs containing an N
    all.ind <- seq(dim(align.umi)[1])
    align.umi <- align.umi[base::setdiff(all.ind, grep("N", align.umi$UMI)), ]
    # tabulate the UMIs
    temp <- as.data.frame(table(align.umi$UMI))
    # keep only UMIs with at least 5 copies but not more than
    align.umi <- align.umi[align.umi$UMI %in% temp[temp[,2] >= min.umi.count & temp[,2] <= max.umi.count,1],]

    # if R2 is within the acceptable length window, then work with the R2 read
    R2.good.len <- subset(align.umi, qwidth.last <= max.R2.len & qwidth.last >= min.R2.mapped)
    R2.umi.plus <- subset(R2.good.len, strand.last == "+")
    R2.umi.minus <- subset(R2.good.len, strand.last == "-")

    # among the leftover reads, check if the R1 is in the acceptable range, and work with that instead. It looks like some of the reads are duplicated; is this intentional?
    R1.good.len <- subset(align.umi, qwidth.first <= max.R1.len & qwidth.first >= min.R1.mapped)
    R1.umi.plus <- subset(R1.good.len, strand.first == "-" & !(readName %in% R2.umi.plus$readName))
    R1.umi.minus <- subset(R1.good.len, strand.first == "+" & !(readName %in% R2.umi.minus$readName))

    #
    unique.umi.plus.R2 <- R2.umi.plus |>
      dplyr::select(seqnames.last, seqnames.first,
             strand.last, strand.first,
             start.last, end.first, UMI) |>
      dplyr::add_count(seqnames.last, seqnames.first,
                strand.last, strand.first,
                start.last, end.first, UMI) |>
      dplyr::distinct() |>
      dplyr::filter(n >= min.read.coverage)

    unique.umi.minus.R2 <- R2.umi.minus %>%
      select(seqnames.last, seqnames.first,
             strand.last, strand.first,
             start.first, end.last, UMI) %>%
      add_count(seqnames.last, seqnames.first,
                strand.last, strand.first,
                start.first, end.last, UMI) %>%
      unique %>%
      filter(n >= min.read.coverage)

    unique.umi.plus.R1 <- R1.umi.plus %>%
      select(seqnames.last, seqnames.first,
             strand.last, strand.first,
             start.first, end.first, UMI) %>%
      add_count(seqnames.last, seqnames.first,
                strand.last, strand.first,
                start.first, end.first, UMI) %>%
      unique %>%
      filter(n >= min.read.coverage)

    unique.umi.minus.R1 <- R1.umi.minus  %>%
      select(seqnames.last, seqnames.first,
             strand.last, strand.first,
             start.first, end.first, UMI) %>%
      add_count(seqnames.last, seqnames.first,
                strand.last, strand.first,
                start.first, end.first, UMI) %>%
      unique %>%
      filter(n >= min.read.coverage)

    plus.cleavage.R2 <-
      unique.umi.plus.R2[, c("seqnames.last", "start.last", "UMI")]
    plus.cleavage.R1 <-
      unique.umi.plus.R1[, c("seqnames.first", "start.first", "UMI")]
    minus.cleavage.R2 <-
      unique.umi.minus.R2[, c("seqnames.last", "end.last", "UMI")]
    minus.cleavage.R1 <-
      unique.umi.minus.R1[, c("seqnames.first", "end.first", "UMI")]
    if (ignoreTagmSite) {
      plus.cleavage.R2 <-  unique(plus.cleavage.R2)
      plus.cleavage.R1 <- unique(plus.cleavage.R1)
      minus.cleavage.R2 <- unique(minus.cleavage.R2)
      minus.cleavage.R1 <- unique(minus.cleavage.R1)
    }

    colnames(plus.cleavage.R1) <- c("seqnames", "start", "UMI")
    colnames(plus.cleavage.R2) <- c("seqnames", "start", "UMI")
    colnames(minus.cleavage.R1) <- c("seqnames", "start", "UMI")
    colnames(minus.cleavage.R2) <- c("seqnames", "start", "UMI")
    plus.cleavage <- rbind(plus.cleavage.R2, plus.cleavage.R1)
    minus.cleavage <- rbind(minus.cleavage.R2, minus.cleavage.R1)
    plus.cleavage <- cbind(plus.cleavage, strand = "+")
    minus.cleavage <- cbind(minus.cleavage, strand = "-")
    colnames(plus.cleavage)[4] <-  "strand"
    colnames(minus.cleavage)[4] <-  "strand"

    unique.umi.both <- rbind(plus.cleavage, minus.cleavage)

    R1.umi.plus <- R1.umi.plus[, c("seqnames.first",
                                   "strand.first",
                                   "start.first", "UMI")]
    R1.umi.minus <- R1.umi.minus[, c("seqnames.first",
                                     "strand.first",
                                     "end.first", "UMI")]
    R2.umi.plus <- R2.umi.plus[, c("seqnames.last",
                                   "strand.last",
                                   "start.last", "UMI")]
    R2.umi.minus <- R2.umi.minus[, c("seqnames.last",
                                     "strand.last",
                                     "end.last", "UMI")]

    colnames(R1.umi.plus) <- c("seqnames", "strand", "start", "UMI")
    colnames(R1.umi.minus) <- c("seqnames", "strand", "start", "UMI")
    colnames(R2.umi.plus) <- c("seqnames", "strand", "start", "UMI")
    colnames(R2.umi.minus) <- c("seqnames", "strand", "start", "UMI")

    # summary before duplicate removal
    R1.umi.plus.summary <- unique(add_count(R1.umi.plus, seqnames,
                                            strand, start, UMI))
    R1.umi.minus.summary <- unique(add_count(R1.umi.minus, seqnames,
                                             strand, start, UMI))
    R2.umi.plus.summary <- unique(add_count(R2.umi.plus, seqnames,
                                            strand, start, UMI))
    R2.umi.minus.summary <- unique(add_count(R2.umi.minus, seqnames,
                                             strand, start, UMI))

    res <- list(cleavage.gr = GRanges(IRanges(
      start=as.numeric(unlist(unique.umi.both[,2])),
      width=1),
      seqnames=unlist(unique.umi.both[,1]),
      strand = unlist(unique.umi.both[,4]),
      total=rep(1, dim(unique.umi.both)[1]),
      umi = unlist(unique.umi.both[,3])),
      unique.umi.plus.R2 = unique.umi.plus.R2,
      unique.umi.minus.R2 = unique.umi.minus.R2,
      unique.umi.plus.R1 = unique.umi.plus.R1,
      unique.umi.minus.R1 = unique.umi.minus.R1,
      align.umi = align.umi,
      umi.count.summary = rbind(R1.umi.plus.summary,
                                R1.umi.minus.summary,
                                R2.umi.plus.summary,
                                R2.umi.minus.summary),
      sequence.depth = length(unique(align$readName))
    )

    res
}


importBAMAlignments <- function(alignment.inputfile,
                                min.mapping.quality = 30L,
                                min.R1.mapped = 20L,
                                min.R2.mapped = 20L,
                                keep.R1only = TRUE,
                                keep.R2only = TRUE,
                                apply.both.min.mapped = FALSE,
                                max.paired.distance = 1000L) {
  # set the loading params; keep only reads whose mapq score exceeds 30
  param <- Rsamtools::ScanBamParam(mapqFilter = min.mapping.quality, what = "flag")
  # load the reads as mates, keeping only those with a mapping (mapq) score >= 30
  gal <- GenomicAlignments::readGAlignmentsList(Rsamtools::BamFile(alignment.inputfile, asMates = TRUE),
                                                param = param, use.names = TRUE)
  # convert to GAlignmentPairs object
  my.pairs <- as(gal, "GAlignmentPairs")
  # filter out chromosome M
  my.pairs <- my.pairs[!seqnames(my.pairs) %in% c("chrM", "chrMT", "M"),]
  # keep only reads whose R1 width or R2 width exceeds 30
  my.pairs <- my.pairs[width(first(my.pairs)) >= min.R1.mapped | width(second(my.pairs)) >= min.R2.mapped]
  # keep only reads where the strand is concordant (i.e., where the R1 and R2 reads map to opposite strands)
  my.pairs <- my.pairs[strand(my.pairs) != "*"]
  # keep only reads where R1 and R2 map to the same chromosome
  my.pairs <- my.pairs[!is.na(seqnames(my.pairs))]
  # the strand of my.pairs is set (by default) to the strand of R1.
  # if R1 maps to +, then R1 is "to the left" of R2. if R2 maps to -, then R1 is "to the right" of R2.
  distance <- ifelse(strand(my.pairs) == "+",
                     start(second(my.pairs)) - end(first(my.pairs)),
                     start(first(my.pairs)) - end(second(my.pairs))
                     )
  mcols(my.pairs)$distance <- distance

  # final step: deal with unpaired reads
  # get the data frame of unpaired reads
  unpaired_df <- unlist(gal[mcols(gal)$mate_status == "unmated"])
  mcols(unpaired_df)$readName <- names(unpaired_df)
  names(unpaired_df) <- NULL
  # determine whether the read is first (i.e., R1)
  first <- Rsamtools::bamFlagTest(mcols(unpaired_df)$flag, "isFirstMateRead")
  # construct the df of R1 unpaired reads
  firstUnpaired_df <- unpaired_df[first & keep.R1only]
  firstUnpaired_df <- firstUnpaired_df[width(firstUnpaired_df) >= min.R1.mapped]
  # construct the df of R2 unpaired reads
  secondUnpaired_df <- unpaired_df[!first & keep.R2only]
  secondUnpaired_df <- secondUnpaired_df[width(secondUnpaired_df) >= min.R2.mapped]
  # combine the R1 and R2 unpaired data frames
  unpairedDF <- rbind_dodge(as.data.frame(firstUnpaired_df),
                            as.data.frame(secondUnpaired_df),
                            ".first", ".last", "readName")
  # set distance to NA -- not needed
  unpairedDF$distance <- NA_integer_
  # convert the my.pairs data frame to pairedDF
  pairedDF <- as.data.frame(my.pairs)
  pairedDF$readName <- rownames(pairedDF)
  rownames(pairedDF) <- NULL

  # combined the paired df with the unpaired df to obtain final data frame
  df <- rbind(pairedDF[colnames(unpairedDF)], unpairedDF)
  df$njunc.first <- df$njunc.last <- NULL

  df
}


rbind_dodge <- function(x, y,
                        xSuffix = deparse(substitute(x)),
                        ySuffix = deparse(substitute(y)),
                        common = character())
{
  naDF <- function(cols) {
    as.data.frame(setNames(as.list(rep(NA, length(cols))), cols))
  }

  dodgeDF <- function(df, suffix) {
    only <- base::setdiff(colnames(df), common)
    dodged <- df[only]
    colnames(dodged) <- paste0(only, suffix)
    dodged
  }

  xDodged <- dodgeDF(x, xSuffix)
  yDodged <- dodgeDF(y, ySuffix)

  rbind(cbind(x[common], xDodged, naDF(colnames(yDodged))),
        cbind(y[common], naDF(colnames(xDodged)), yDodged))
}
