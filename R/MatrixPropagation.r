#' Propagation of the phenotypic score within a multilayer network with
#' physical and functional interactions
#'
#' @description
#' Propagates the phenotypic prior score across two independent gene interaction
#' networks (physical and functional) using an iterative random-walk algorithm.
#' At each iteration, each gene's score is updated as a weighted combination of
#' its neighbours' scores (\code{alpha}) and its original prior score
#' (\code{1 - alpha}). Convergence is assessed by rank-order stability over a
#' sliding window of recent iterations. Returns a ranked gene table that is used
#' as \code{GlobalPhenotypicScore} in \code{priorBestVariant()}.
#'
#' @param Y Numeric vector with the phenotypic prior score for each gene, as
#'   returned by \code{proteinScore()$scores}.
#' @param alpha Numeric [0, 1]. Propagation weight: balances network diffusion
#'   (\code{alpha}) vs. retention of the original prior score (\code{1 - alpha})
#'   at each iteration. Higher values allow more diffusion. Default \code{0.2}.
#' @param max_iter Integer. Maximum number of iterations before stopping,
#'   regardless of convergence. Default \code{10000}.
#' @param window Integer. Number of consecutive iterations with identical rank
#'   order required to declare convergence. Default \code{20}.
#' @param verbose Logical. If \code{TRUE} (default), prints the number of
#'   iterations at convergence for each network.
#'
#' @return A \code{data.frame} with one row per gene and six columns:
#'   \code{Symbol_Funct} and \code{geneID_Funct} (gene symbol and Ensembl ID,
#'   sorted by functional score), \code{PriorFunct} (propagated functional
#'   network score), \code{Symbol_Phys} and \code{geneID_Phys} (sorted by
#'   physical score), and \code{PriorPhys} (propagated physical network score).
#'   The two networks are ranked independently. Pass this object directly as
#'   \code{GlobalPhenotypicScore} to \code{priorBestVariant()}.
#'
#' @export
#' @importFrom Matrix Matrix rowSums
#'
#' @examples
#' HPOpatient <- c("HP:0004481", "HP:0002376", "HP:0001257", "HP:0001250",
#'                 "HP:0000238", "HP:0002922", "HP:0000365")
#' result              <- proteinScore(HPOpatient)
#' GlobalPhenotypicScore <- MatrixPropagation(result$scores, alpha = 0.2)
MatrixPropagation <- function(Y, alpha = 0.2, max_iter = 10000,
                              window = 20, verbose = TRUE) {
  
  # -- Input validation -------------------------------------------------------
  if (!is.numeric(Y))            stop("Y must be a numeric vector.")
  if (alpha < 0 || alpha > 1)   stop("alpha must be between 0 and 1.")
  if (length(Y) != nrow(total_unique))
    warning("Length of Y (", length(Y), ") differs from nrow(total_unique) (",
            nrow(total_unique), "). Check inputs.", call. = FALSE)
  
  # -- Internal helper: iterative propagation ---------------------------------
  propagate <- function(F_init, pos, norm_matrix, network_name) {
    F       <- F_init
    F2      <- F[pos]
    memo    <- F2
    memoL   <- rep(NA_integer_, window)   # circular buffer for convergence
    iter    <- 0L
    
    repeat {
      iter  <- iter + 1L
      WF    <- norm_matrix %*% F2
      F2    <- (alpha * WF) + (1 - alpha) * F[pos]
      
      # Convergence: compare rank order of current vs previous iteration
      l     <- sum(names(sort(memo)) == names(sort(F2)))   # rank-stable count
      # fallback if F2 has no names (use positional comparison)
      l     <- sum(order(as.numeric(memo)) == order(as.numeric(F2)))
      
      memoL <- c(memoL[-1], l)   # slide the window
      
      converged <- !anyNA(memoL) && all(memoL == l)
      if (converged || iter >= max_iter) {
        if (verbose) message(
          "[MatrixPropagation] ", network_name,
          " converged after ", iter, " iteration(s)."
        )
        if (iter >= max_iter && !converged)
          warning("[MatrixPropagation] ", network_name,
                  " reached max_iter (", max_iter, ") without full convergence.",
                  call. = FALSE)
        break
      }
      memo <- F2
    }
    
    F[pos] <- F2
    F
  }
  
  # -- Physical network -------------------------------------------------------
  pos_phys  <- match(total_unique_Conn_Physical, total_unique[, 2])
  if (anyNA(pos_phys))
    warning(sum(is.na(pos_phys)),
            " gene(s) in total_unique_Conn_Physical not found in total_unique.",
            call. = FALSE)
  
  F_phys    <- propagate(Y, pos_phys, normPhysical, "Physical")
  ord_phys  <- order(F_phys, decreasing = TRUE)
  result_phys <- data.frame(
    Symbol_Phys  = total_unique[ord_phys, 1],
    geneID_Phys  = total_unique[ord_phys, 2],
    PriorPhys    = as.numeric(F_phys[ord_phys]),
    stringsAsFactors = FALSE
  )
  
  # -- Functional network -----------------------------------------------------
  pos_func  <- match(total_unique_Conn_Func, total_unique[, 2])
  if (anyNA(pos_func))
    warning(sum(is.na(pos_func)),
            " gene(s) in total_unique_Conn_Func not found in total_unique.",
            call. = FALSE)
  
  F_func    <- propagate(Y, pos_func, normFunc, "Functional")
  ord_func  <- order(F_func, decreasing = TRUE)
  result_func <- data.frame(
    Symbol  = total_unique[ord_func, 1],
    geneID  = total_unique[ord_func, 2],
    PriorFunct = as.numeric(F_func[ord_func]),
    stringsAsFactors = FALSE
  )
  
  # -- Combine & return -------------------------------------------------------
  # Each network is sorted independently by its own score
  ClinPriorScore <- data.frame(
    Symbol_Funct = total_unique[ord_func, 1],
    geneID_Funct = total_unique[ord_func, 2],
    PriorFunct   = as.numeric(F_func[ord_func]),
    Symbol_Phys  = total_unique[ord_phys, 1],
    geneID_Phys  = total_unique[ord_phys, 2],
    PriorPhys    = as.numeric(F_phys[ord_phys]),
    stringsAsFactors = FALSE
  )
  rownames(ClinPriorScore) <- NULL
  
  gc()
  return(ClinPriorScore)
}