# BioPlex interaction networks (affinity purification-mass spectrometry),
# as published by BioPlex (already filtered by the BioPlex pipeline):
#   293T:   BioPlex 3.0, HEK293T cells   (BioPlex_293T_Network_10K_Dec_2019.tsv)
#   HCT116: BioPlex HCT116 1.0           (BioPlex_HCT116_Network_5.5K_Dec_2019.tsv)
# Columns: GeneA GeneB (Entrez) UniprotA UniprotB SymbolA SymbolB pW pNI pInt.
# The original ClinPrior physical network includes BioPlex 2.0; BioPlex 3.0
# (HEK293T) contains it.

bioplexUrls <- c(
  "293T"   = "https://bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv",
  "HCT116" = "https://bioplex.hms.harvard.edu/data/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv"
)

# files: named vector of local files (names = cell lines).
# min_pint: additional minimum interaction probability (0 = keep all).
readBioPlex <- function(files, min_pint = 0) {
  out <- lapply(names(files), function(cl) {
    b <- utils::read.delim(files[[cl]], sep = "\t", header = TRUE, quote = "\"",
                           stringsAsFactors = FALSE)
    if(ncol(b) < 9) stop("Unexpected BioPlex format in ", files[[cl]])
    names(b)[1:9] <- c("GeneA", "GeneB", "UniprotA", "UniprotB",
                       "SymbolA", "SymbolB", "pW", "pNI", "pInt")
    a    <- suppressWarnings(as.integer(b$GeneA))
    z    <- suppressWarnings(as.integer(b$GeneB))
    pint <- suppressWarnings(as.numeric(b$pInt))
    keep <- !is.na(a) & !is.na(z) & a != z & !is.na(pint) & pint >= min_pint
    list(pairs = data.frame(a = pmin(a, z), b = pmax(a, z))[keep, ],
         symbols = data.frame(entrez = c(a[keep], z[keep]),
                              symbol = c(b$SymbolA[keep], b$SymbolB[keep]),
                              stringsAsFactors = FALSE),
         n = sum(keep))
  })
  pairs <- unique(do.call(rbind, lapply(out, `[[`, "pairs")))
  sym   <- do.call(rbind, lapply(out, `[[`, "symbols"))
  attr(pairs, "symbols") <- sym[!duplicated(sym$entrez), ]
  attr(pairs, "release") <- paste0(names(files), " (", vapply(out, `[[`, 1L, "n"), " interactions)",
                                   collapse = ", ")
  pairs
}