#!/usr/bin/env Rscript

# Load required library
library('BSgenomeForge')

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments are provided
if (length(args) < 6) {
  stop("Usage: forge_bsgenome.R <path> <genus> <species> <assembly> <mt> <library_dir>")
}


# Assign variables from command-line arguments
path_2bit <- args[1]
genus <- args[2]
species <- args[3]
assembly <- args[4]
mt <- args[5]
library_dir <- args[6]



# Define package destination
pkg_name <- paste0("BSgenome.", paste0(toupper(substr(genus, 1, 1)), species), ".NA.", assembly)

# Generate the BSgenome data package
forgeBSgenomeDataPkgFromTwobitFile(path_2bit, paste(genus, species), 'NA', assembly, 
    'NA <nully@null.com>',
    pkg_author = 'NA',
    pkg_version = "1.0.0",
    pkg_license = "Artistic-2.0",
    circ_seqs = c(mt),
    destdir = "." )


# Install the generated package
install.packages(pkg_name, repos = NULL, type = 'source', lib = library_dir)

