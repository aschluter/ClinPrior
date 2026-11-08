# Add BioGRID multi-validated (MV) physical interactions to the physical
# network (normPhysical).
#
# normPhysical is D^-1 (A + I): A holds integer, symmetric edge weights
# (number of times the edge appears in the sources, usually 2) and I a self
# loop of weight 1; each row is divided by its sum. A is recovered exactly
# from the diagonal (weight 1), BioGRID edges not already present are added
# with weight `weight`, and the matrix is normalised again.

# Physical systems that are not protein-protein interactions
.biogrid_not_ppi <- c("Co-localization", "Affinity Capture-RNA", "Protein-RNA")
# Binary (direct) protein-protein systems
.biogrid_binary  <- c("Two-hybrid", "PCA", "Reconstituted Complex", "Co-crystal Structure",
                      "FRET", "Far Western", "Protein-peptide", "Biochemical Activity")

# Human, physical, non-self interactions from a BioGRID tab3 file (or .zip).
# systems: "ppi" (all physical except .biogrid_not_ppi), "binary", or a
# character vector of BioGRID experimental systems.
readBioGRID <- function(file, systems = "ppi") {
  if(grepl("\\.zip$", file)) {
    inner <- utils::unzip(file, list = TRUE)$Name
    inner <- inner[grepl("tab3\\.txt$", inner)][1]
    con   <- unz(file, inner)
    release <- sub("\\.tab3\\.txt$", "", inner)
  } else {
    con <- file
    release <- sub("\\.tab3\\.txt$", "", basename(file))
  }
  b <- utils::read.delim(con, sep = "\t", header = TRUE, quote = "", comment.char = "",
                         check.names = FALSE, stringsAsFactors = FALSE, na.strings = "-")
  col <- function(n) b[[n]]
  keep <- col("Organism ID Interactor A") %in% 9606 & col("Organism ID Interactor B") %in% 9606 &
    col("Experimental System Type") == "physical"
  sys <- col("Experimental System")
  keep <- keep & if(identical(systems, "ppi")) !sys %in% .biogrid_not_ppi else
    if(identical(systems, "binary")) sys %in% .biogrid_binary else sys %in% systems
  a <- suppressWarnings(as.integer(col("Entrez Gene Interactor A")))
  z <- suppressWarnings(as.integer(col("Entrez Gene Interactor B")))
  keep <- keep & !is.na(a) & !is.na(z) & a != z
  out <- data.frame(a = pmin(a, z), b = pmax(a, z))[keep, ]
  sym <- unique(data.frame(entrez = c(a, z)[c(keep, keep)],
                           symbol = c(col("Official Symbol Interactor A"),
                                      col("Official Symbol Interactor B"))[c(keep, keep)]))
  out <- unique(out)
  attr(out, "symbols") <- sym[!duplicated(sym$entrez), ]
  attr(out, "release") <- release
  out
}

# W: row-normalised network; genes: Entrez IDs of its rows/columns;
# pairs: data.frame(a, b) of Entrez IDs.
# Genes of the original network are never removed. New edges are not added
# to genes that would end up with more than max_degree neighbours (they keep
# their original edges); new genes in that situation are not added.
addEdgesToNetwork <- function(W, genes, pairs, weight = 2, max_degree = 1000) {
  # Recover integer weights: each row divided by its diagonal (weight 1)
  Aw <- Matrix::Diagonal(x = 1 / Matrix::diag(W)) %*% W
  Aw@x <- round(Aw@x)
  Aw <- methods::as(Aw, "CsparseMatrix")
  n0 <- length(genes)
  
  # Candidate edges: pairs not already in the network
  cand_genes <- as.integer(c(genes, setdiff(unique(c(pairs$a, pairs$b)), genes)))
  i <- match(pairs$a, cand_genes); j <- match(pairs$b, cand_genes)
  inside  <- i <= n0 & j <= n0
  present <- rep(FALSE, length(i))
  present[inside] <- Aw[cbind(i[inside], j[inside])] != 0
  i <- i[!present]; j <- j[!present]
  
  # Degree each gene would have with all candidate edges
  deg0  <- c(Matrix::rowSums(Aw != 0) - 1, rep(0, length(cand_genes) - n0))
  deg   <- deg0 + tabulate(c(i, j), nbins = length(cand_genes))
  hub   <- deg > max_degree
  keep  <- !hub[i] & !hub[j]
  i <- i[keep]; j <- j[keep]
  
  # New genes that keep at least one edge
  used      <- sort(unique(c(i, j)[c(i, j) > n0]))
  all_genes <- c(genes, cand_genes[used])
  n         <- length(all_genes)
  remap     <- integer(length(cand_genes)); remap[seq_len(n0)] <- seq_len(n0)
  remap[used] <- n0 + seq_along(used)
  if(n > n0) {
    Aw <- rbind(cbind(Aw, Matrix::Matrix(0, n0, n - n0, sparse = TRUE)),
                cbind(Matrix::Matrix(0, n - n0, n0, sparse = TRUE),
                      Matrix::Diagonal(n - n0)))
  }
  if(length(i) > 0) {
    E <- Matrix::sparseMatrix(i = c(remap[i], remap[j]), j = c(remap[j], remap[i]),
                              x = weight, dims = c(n, n))
    Aw <- Aw + E
  }
  
  Wn <- Matrix::Diagonal(x = 1 / Matrix::rowSums(Aw)) %*% Aw
  Wn <- methods::as(Wn, "CsparseMatrix")
  dimnames(Wn) <- list(as.character(all_genes), as.character(all_genes))
  list(W = Wn, genes = all_genes,
       added_edges = length(i), present_edges = sum(present),
       new_genes = cand_genes[used],
       capped = cand_genes[hub])
 }