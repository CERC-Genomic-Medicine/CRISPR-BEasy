#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(jsonlite)
  library(BSgenomeForge)
  library(Biostrings)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    "Usage: Rscript forge_bsgenome.R <species> <Assembly> <FastaFile> <GFF file> <Source>\n"
  )
}


get_circular_seqs <- function(gff3_path, twobit_path) {
  # import GFF
  gff <- rtracklayer::import(gff3_path)

  # logical tests (NA → dropped by which)
  idx <- which(
    mcols(gff)$type        == "region" &
    mcols(gff)$Is_circular == "true"
  )

  # pull unique seqnames of circular regions
  circ_ids <- unique(
    as.character( GenomicRanges::seqnames(gff)[idx] )
  )

  # open 2bit and get its seqnames
  tb      <- rtracklayer::TwoBitFile(twobit_path)
  tb_ids  <- as.character( GenomeInfoDb::seqnames(seqinfo(tb)) )

  # keep exact matches and mitochondrial aliases
  mito_aliases <- tolower(c("MT","mt","chrMT","M","chrM"))
  keep_circ    <- tb_ids[ tb_ids %in% circ_ids ]
  keep_mito    <- tb_ids[ tolower(tb_ids) %in% mito_aliases ]

  sort(unique(c(keep_circ, keep_mito)))
}


SpeciesUC <- args[1]  # e.g. "Homo_sapiens"
Assembly  <- gsub("\\.", "_",args[2])  # e.g. "GRCh38"
path_2bit <- normalizePath(args[3], mustWork = TRUE)
GFFfile <- normalizePath(args[4], mustWork = TRUE)
sourced <- args[5]

SpeciesUC <- paste0(substr(SpeciesUC, 1, 1), sub(".*_", "", SpeciesUC))
circ_seqs <- get_circular_seqs(GFFfile,path_2bit)

forgeBSgenomeDataPkgFromTwobitFile(path_2bit, SpeciesUC, sourced, Assembly,
    'NA <nully@null.com>',
    pkg_author = 'NA',
    pkg_version = "1.0.0",
    pkg_license = "Artistic-2.0",
    circ_seqs = circ_seqs,
    destdir = "." )
