#' Compute phenotypic gene scores from patient HPO terms
#'
#' @description
#' Computes a phenotypic prior score for each gene based on the patient's
#' clinical description in HPO terms. The function validates and filters the
#' input HPO terms, builds a subgraph of phenotypically related terms, groups
#' them by cluster, and applies a hypergeometric test to quantify the
#' gene-phenotype match. The resulting scores are normalised and transformed
#' via a sigmoid function. The output is designed to be passed directly to
#' \code{MatrixPropagation()} as \code{result$scores}.
#'
#' @param HPOpatient Character vector with the patient's HPO term identifiers
#'   (e.g. \code{c("HP:0001250", "HP:0001249")}). Duplicates are removed
#'   automatically.
#' @param verbose Logical. If \code{TRUE} (default), prints diagnostic messages
#'   about unrecognized, removed, and used HPO terms.
#'
#' @return A list with four elements:
#'   \describe{
#'     \item{scores}{Numeric vector of length equal to \code{nrow(total_unique)}
#'       with the sigmoid-normalised phenotypic score for each gene. Pass this
#'       to \code{MatrixPropagation()} as the \code{Y} argument.}
#'     \item{unrecognized_HPOs}{Character vector of HPO terms not found in the
#'       internal HPO adjacency matrix (\code{HPOadj}).}
#'     \item{removed_HPOs}{Character vector of HPO terms excluded as irrelevant
#'       based on the internal \code{treureHPO} filter list.}
#'     \item{used_HPOs}{Character vector of HPO terms actually used for scoring,
#'       after removing unrecognized and irrelevant terms.}
#'   }
#'
#' @export
#' @importFrom igraph as_undirected ego induced_subgraph cluster_edge_betweenness
#' @importFrom Matrix Matrix rowSums
#' @importFrom stats phyper
#'
#' @examples
#' HPOpatient <- c("HP:0004481", "HP:0002376", "HP:0001257", "HP:0001250",
#'                 "HP:0000238", "HP:0002922", "HP:0000365")
#' result <- proteinScore(HPOpatient)
#' result$scores           # numeric vector of gene scores
#' result$unrecognized_HPOs  # HPO terms not found in the database
#' result$used_HPOs          # HPO terms effectively used for scoring
proteinScore <- function(HPOpatient, verbose = TRUE) {
  
  # -- 1. Input validation & deduplication -----------------------------------
  HPOpatient <- unique(as.character(HPOpatient))
  
  if (length(HPOpatient) == 0) stop("HPOpatient is empty.")
  
  # -- 2. Track unrecognized HPOs --------------------------------------------
  # An HPO is "recognized" if it appears in the adjacency matrix columns
  known_HPOs       <- colnames(HPOadj)
  unrecognized     <- HPOpatient[!HPOpatient %in% known_HPOs]
  HPOpatient_valid <- HPOpatient[HPOpatient %in% known_HPOs]
  
  if (length(unrecognized) > 0) {
    msg <- paste0(
      "[proteinScore] ", length(unrecognized),
      " HPO term(s) not recognized (not in HPOadj):\n  ",
      paste(unrecognized, collapse = ", ")
    )
    if (verbose) message(msg)
    warning(msg, call. = FALSE)
  }
  
  if (length(HPOpatient_valid) == 0) {
    stop("No valid HPO terms remain after removing unrecognized terms.")
  }
  
  # -- 3. Remove irrelevant HPOs ---------------------------------------------
  treureHPO_pattern <- paste(treureHPO, collapse = "|")
  
  # From the gene-HPO table
  HPO2genes_filt  <- HPO2genes
  posTreure_genes <- grep(treureHPO_pattern, HPO2genes_filt[, 2])
  if (length(posTreure_genes) > 0) HPO2genes_filt <- HPO2genes_filt[-posTreure_genes, ]
  
  # From the patient HPO list
  posTreure_patient <- grep(treureHPO_pattern, HPOpatient_valid)
  removed_HPOs      <- character(0)
  if (length(posTreure_patient) > 0) {
    removed_HPOs     <- HPOpatient_valid[posTreure_patient]
    HPOpatient_valid <- HPOpatient_valid[-posTreure_patient]
    if (verbose) message(
      "[proteinScore] ", length(removed_HPOs),
      " HPO term(s) removed as irrelevant:\n  ",
      paste(removed_HPOs, collapse = ", ")
    )
  }
  
  if (length(HPOpatient_valid) == 0) {
    stop("No valid HPO terms remain after removing irrelevant terms.")
  }
  
  used_HPOs <- HPOpatient_valid
  if (verbose) message(
    "[proteinScore] ", length(used_HPOs),
    " HPO term(s) used for scoring:\n  ",
    paste(used_HPOs, collapse = ", ")
  )
  
  # -- 4. Build subgraph & cluster -------------------------------------------
  g1_undir <- igraph::as_undirected(g1)                          
  
  HPOorig_expanded <- unique(unlist(
    sapply(used_HPOs, function(x) {
      ego_result <- igraph::ego(g1_undir, order = 1, nodes = x)
      if (length(ego_result) == 0 || is.null(ego_result[[1]])) return(NULL)
      rownames(as.matrix(ego_result[[1]]))
    })
  ))
  
  g.sub <- igraph::induced_subgraph(graph = g1_undir,   
                                    vids  = HPOorig_expanded)
  res                  <- cluster_edge_betweenness(g.sub)
  HPOorig_expanded_mat <- cbind(res$names, res$membership)
  
  # Align groups to the (valid, non-removed) patient HPOs
  HPOorigGroups <- HPOorig_expanded_mat[
    match(used_HPOs, HPOorig_expanded_mat[, 1]), , drop = FALSE
  ]
  
  # -- 5. Accumulate scores --------------------------------------------------
  genes    <- unique(HPO2genes_filt[, 1])
  posGenes <- match(genes, rownames(HPOadj))
  
  acumulat <- Matrix(matrix(0, nrow = length(posGenes), ncol = length(used_HPOs)))
  
  for (z in seq_along(used_HPOs)) {
    pos            <- match(used_HPOs[z], colnames(HPOadj))
    HPOPatientItem <- HPOdistance[pos, ]
    acumulat[, z]  <- HPOadj[posGenes, ] %*% HPOPatientItem
  }
  
  # -- 6. Merge HPOs belonging to the same cluster group --------------------
  groups <- HPOorigGroups[, 2]
  if (length(groups) != length(unique(groups))) {
    dupli        <- unique(groups[duplicated(groups)])
    cols_to_drop <- c()
    for (i in seq_along(dupli)) {
      pos                <- which(groups == dupli[i])
      acumulat[, pos[1]] <- Matrix::rowSums(acumulat[, pos, drop = FALSE])
      cols_to_drop       <- c(cols_to_drop, pos[-1])
    }
    acumulat <- acumulat[, -cols_to_drop, drop = FALSE]
  }
  
  # -- 7. Binarise & count matched groups -----------------------------------
  acumulat[is.na(acumulat)] <- 0
  acumulat                  <- Matrix(acumulat / acumulat)
  acumulat[is.na(acumulat)] <- 0
  
  HPOmatch_quant <- Matrix::rowSums(acumulat)
  
  # -- 8. Hypergeometric test ------------------------------------------------
  q <- HPOmatch_quant - 1
  m <- length(unique(HPOorigGroups[, 2]))
  n <- length(igraph::V(g1)$name) - m
  k <- as.matrix(HPOqueryGene) * 10
  
  stats <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  gc()
  
  stats[stats == -Inf] <- 1
  min_nonzero          <- min(stats[stats != 0], na.rm = TRUE)
  stats[stats == 0]    <- 10^(log10(min_nonzero) - 1)
  
  # -- 9. Score normalisation & sigmoid -------------------------------------
  D      <- abs(log10(abs(stats)))
  Dred   <- as.numeric(D)
  
  pos    <- match(genes, total_unique[, 2])
  D_full <- numeric(nrow(total_unique))
  D_full[pos] <- Dred
  
  DNormed <- (D_full - min(D_full, na.rm = TRUE)) /
    (max(D_full, na.rm = TRUE) - min(D_full, na.rm = TRUE))
  
  Y <- as.numeric(1 / (1 + exp((DNormed * (-12)) + log(9999))))
  
  # -- 10. Return ------------------------------------------------------------
  return(list(
    scores            = Y,
    unrecognized_HPOs = unrecognized,
    removed_HPOs      = removed_HPOs,
    used_HPOs         = used_HPOs
  ))
}