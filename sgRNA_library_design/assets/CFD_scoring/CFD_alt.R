#!/usr/bin/env Rscript

PAM_to_DICT <- function(pattern) {
  # Convert pattern to uppercase.
  pattern <- toupper(pattern)
  
  # Define mapping from ambiguous codes to allowed nucleotides.
  mapping <- list(
  A = c("A"),
  C = c("C"),
  G = c("G"),
  T = c("T"),
  U = c("U"),
  R = c("A", "G"),             # purine
  Y = c("C", "T"),             # pyrimidine
  S = c("G", "C"),
  W = c("A", "T"),
  K = c("G", "T"),
  M = c("A", "C"),
  B = c("C", "G", "T"),
  D = c("A", "G", "T"),
  H = c("A", "C", "T"),
  V = c("A", "C", "G"),
  N = c("A", "C", "G", "T")
  )
  
  # Split the pattern into individual characters.
  pattern_chars <- unlist(strsplit(pattern, split = ""))
  
  # Check for unrecognized nucleotide codes.
  invalid_codes <- pattern_chars[!pattern_chars %in% names(mapping)]
  if (length(invalid_codes) > 0) {
    stop(paste("Unrecognized nucleotide code in pattern:",
               paste(unique(invalid_codes), collapse = ", ")))
  }
  
  # Create a list of allowed bases for each position in the pattern.
  allowed_bases <- lapply(pattern_chars, function(x) mapping[[x]])
  
  # Generate all concrete sequences that match the input pattern.
  # Use expand.grid to compute the Cartesian product of allowed bases.
  valid_df <- do.call(expand.grid, c(allowed_bases, stringsAsFactors = FALSE))
  valid_seq <- apply(valid_df, 1, paste0, collapse = "")
  
  # Generate all possible combinations of A, T, C, G for sequences of this length.
  n <- nchar(pattern)
  all_bases <- rep(list(c("A", "T", "C", "G")), n)
  all_df <- do.call(expand.grid, c(all_bases, stringsAsFactors = FALSE))
  all_seq <- apply(all_df, 1, paste0, collapse = "")
  
  # Create a named numeric vector:
  # 1.0 for sequences that fulfill the pattern, 0.0 otherwise.
  result <- ifelse(all_seq %in% valid_seq, 1.0, 0.0)
  names(result) <- all_seq
  
  return(result)
}

getCFDScores_ALT_cas9 <- function(spacers,
                         protospacers,
                         pams,
                         pam_weights,
                         Mismatch_dict)
    {
    if (length(spacers)==1){
        spacers <- rep(spacers, length(protospacers))
    } else {
        if (length(spacers) != length(protospacers)){
            stop("Input vectors 'spacers' and 'protospacers' must",
                 " have the same length.")
        }
    }

    if (unique(nchar(protospacers))!=20){
        stop("Protospacer sequences must have length 20nt.")
    }
    if (unique(nchar(spacers))!=20){
        stop("Spacer sequences must have length 20nt.")
    }

    spacers.wt  <- spacers
    spacers.off <- protospacers
    spacers.wt   <- DNAStringSet(spacers.wt)
    spacers.off  <- DNAStringSet(spacers.off)
    x <- as.matrix(spacers.wt)
    y <- as.matrix(spacers.off)
    wh <- x!=y
    spacer.len <- unique(width(spacers.wt))
    MM <- matrix(rep(seq_len(spacer.len), each=nrow(x)),
                 nrow=nrow(x), ncol=ncol(x))
    MM[wh] <- paste0(x[wh], y[wh], MM[wh])
    mm <- lapply(seq_len(nrow(MM)), function(i){
        x <- MM[i,]
        xx <- x[grepl("[ACGT]+",x)]
        if (length(xx)==0){
            xx <- NA
        }
        return(xx)
    })

    mm.scores <- vapply(mm,
                        function(x) prod(Mismatch_dict[x]),
                        FUN.VALUE=0)
    mm.scores[is.na(mm.scores)] <- 1
    pam.scores <- pam_weights[pams]
    cfd.scores <- as.numeric(mm.scores*pam.scores)
    data.frame(spacer=spacers,
               protospacer=protospacers,
               score=cfd.scores)
    return(cfd.scores)
}
