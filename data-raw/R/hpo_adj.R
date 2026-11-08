# Gene x HPO matrices used by proteinScore().
#
# HPOadj:       genes (rows) x HPO terms of g1 (columns), 1 if annotated.
#               Rows: unique genes of HPO2genes after removing treureHPO, in
#               table order (proteinScore() relies on this order).
# HPOdistance:  HPO x HPO, 1 if the terms are the same or directly linked
#               in g1 (shortest undirected distance < 2). Built from the
#               sparse adjacency matrix, which gives the same matrix as
#               igraph::distances() without the dense HPO x HPO matrix.
# HPOqueryGene: number of annotated HPO terms per gene (rowSums of HPOadj).
buildHPOadj <- function(HPO2genes, g1, treureHPO) {
  HPO_terms <- igraph::V(g1)$name
  
  posTreure <- grep(paste(treureHPO, collapse = "|"), HPO2genes[, 2])
  if(length(posTreure) > 0) HPO2genes <- HPO2genes[-posTreure, ]
  
  genes <- unique(HPO2genes[, 1])
  valid <- !is.na(match(HPO2genes[, 2], HPO_terms))
  h     <- HPO2genes[valid, ]
  
  HPOadj <- Matrix::sparseMatrix(
    i = match(h[, 1], genes), j = match(h[, 2], HPO_terms), x = 1,
    dims = c(length(genes), length(HPO_terms)),
    dimnames = list(as.character(genes), HPO_terms))
  HPOadj@x[] <- 1   # duplicated pairs count once
  
  A <- igraph::as_adjacency_matrix(g1, sparse = TRUE)
  HPOdistance <- Matrix::forceSymmetric(
    Matrix::drop0((A + Matrix::t(A) + Matrix::Diagonal(nrow(A))) > 0) * 1)
  dimnames(HPOdistance) <- list(HPO_terms, HPO_terms)
  
  HPOqueryGene <- HPOadj %*% rep(1, ncol(HPOadj))
  
  list(HPOadj = HPOadj, HPOdistance = HPOdistance, HPOqueryGene = HPOqueryGene)
}