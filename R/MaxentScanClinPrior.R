#' Evaluate splice site strength changes using MaxEntScan (worker function)
#'
#' @description
#' Internal worker function executed in each parallel process by
#' \code{MaxentScanClinPrior()}. For each variant, extracts the genomic sequence
#' context and computes MaxEntScan scores for donor (5'ss) and acceptor (3'ss)
#' splice sites, comparing reference and alternative alleles. Handles both
#' forward and reverse strand genes.
#'
#' @param idx Integer vector. Row indices of \code{matrixVariants} to process
#'   in this chunk.
#' @param matrixVariants A \code{data.frame} with columns \code{CHROM},
#'   \code{POS}, \code{REF}, \code{ALT}, and \code{genesList}.
#' @param assembly Character. Genome assembly version. Either \code{"assembly38"}
#'   (hg38) or \code{"assembly37"} (hg19).
#' @param Threshold Numeric. Minimum fractional score change to consider a
#'   splice site as affected. Default \code{0.2}.
#' @param maxentpyPath Character. Path to the \code{maxentpy} Python module,
#'   as returned by \code{system.file("extdata", package = "ClinPrior")}.
#'
#' @return A matrix with the input variant columns plus 8 additional columns:
#'   \code{acceptor_loss}, \code{acceptor_loss_pos}, \code{acceptor_gain},
#'   \code{acceptor_gain_pos}, \code{donor_loss}, \code{donor_loss_pos},
#'   \code{donor_gain}, \code{donor_gain_pos}.
#'
#' @importFrom Biostrings DNAString reverseComplement
#' @importFrom BSgenome getSeq
#'
#' @keywords internal
.MaxentScanWorker <- function(idx, matrixVariants, assembly, Threshold, maxentpyPath) {
  
  matrixVariants_chunk <- matrixVariants[idx, ]
  
  # Carreguem en l'ordre correcte - reticulate SEMPRE l'ultim
  # per evitar conflictes amb BiocIO/rtracklayer
  requireNamespace("Biostrings",    quietly = TRUE)
  requireNamespace("BSgenome",      quietly = TRUE)
  requireNamespace("reticulate",    quietly = TRUE)
  
  BSgenome_assembly <- if (assembly == "assembly37") {
    requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
    BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else {
    requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
    BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }
  
  maxentpy <- reticulate::import_from_path("maxentpy", path = maxentpyPath, convert = TRUE)
  maxentpy$maxent$dir_path <- maxentpyPath
  
  reticulate::py_run_string("
import maxentpy.maxent as maxent

def score3_batch(seqs):
    results = []
    for s in seqs:
        try:
            results.append(maxent.score3(str(s)))
        except Exception:
            results.append(0.0)
    return results

def score5_batch(seqs):
    results = []
    for s in seqs:
        try:
            results.append(maxent.score5(str(s)))
        except Exception:
            results.append(0.0)
    return results
")
  
  score3_batch <- function(seqs) as.numeric(unlist(reticulate::py$score3_batch(as.list(seqs))))
  score5_batch <- function(seqs) as.numeric(unlist(reticulate::py$score5_batch(as.list(seqs))))
  
  getStrand <- function(gene) gene2strand[match(gene, gene2strand[, 1]), 2]
  
  # -------------------------------------------------------
  # Ref2Alt3 (acceptor 3'ss): retorna score i posicio
  # -------------------------------------------------------
  Ref2Alt3 <- function(seq1, seq2, change1, change2) {
    m    <- gregexpr("[A|C|G|T]{18}AG[A|C|G|T]{3}", as.character(seq1), perl = TRUE)
    seqs <- regmatches(as.character(seq1), m)[[1]]
    
    if (length(seqs) == 0) return(list(list(score = 0, pos = NA)))
    
    lapply(seqs, function(seqOne) {
      score <- 0
      pos   <- NA
      l     <- gregexpr(seqOne, as.character(seq1), perl = TRUE)[[1]][1]
      
      if (l > 3 && l < (27 + nchar(change1))) {
        score <- score3_batch(seqOne)
        if (score >= 3) {
          altSeq <- as.character(seq2[l:(l + 22)])
          score2 <- score3_batch(altSeq)
          if (score2 < score) {
            score <- abs(score - score2) / score
            pos   <- l
          } else {
            score <- 0
          }
          if (score < Threshold) { score <- 0; pos <- NA }
        } else {
          score <- 0
        }
      }
      list(score = score, pos = pos)
    })
  }
  
  # -------------------------------------------------------
  # Ref2Alt5 (donor 5'ss): retorna score i posicio
  # -------------------------------------------------------
  Ref2Alt5 <- function(seq1, seq2, change1, change2) {
    m    <- gregexpr("[A|C|G|T]{3}GT[A|C|G|T]{4}", as.character(seq1), perl = TRUE)
    seqs <- regmatches(as.character(seq1), m)[[1]]
    
    if (length(seqs) == 0) return(list(list(score = 0, pos = NA)))
    
    lapply(seqs, function(seqOne) {
      score <- 0
      pos   <- NA
      l     <- gregexpr(seqOne, as.character(seq1), perl = TRUE)[[1]][1]
      
      if (l > 18 && l < (27 + nchar(change1))) {
        score <- score5_batch(seqOne)
        if (score >= 3) {
          altSeq <- as.character(seq2[l:(l + 8)])
          score2 <- score5_batch(altSeq)
          if (score2 < score) {
            score <- abs(score - score2) / score
            pos   <- l
          } else {
            score <- 0
          }
          if (score < Threshold) { score <- 0; pos <- NA }
        } else {
          score <- 0
        }
      }
      list(score = score, pos = pos)
    })
  }
  
  # -------------------------------------------------------
  # Helper: millor score i posicio
  # -------------------------------------------------------
  getBestScoreAndPos <- function(results) {
    scores <- sapply(results, function(x) x$score)
    best_i <- order(unlist(scores), decreasing = TRUE)[1]
    list(
      score = unlist(scores)[best_i],
      pos   = results[[best_i]]$pos
    )
  }
  
  relToGenomic <- function(pos_rel, variant_pos, type) {
    if (is.na(pos_rel)) return(NA)
    offset <- if (type == "acceptor") 18 else 3
    as.numeric(variant_pos) - 26 + pos_rel + offset
  }
  
  relToGenomicMinus <- function(pos_rel, variant_pos, type) {
    if (is.na(pos_rel)) return(NA)
    offset <- if (type == "acceptor") 18 else 3
    as.numeric(variant_pos) + 26 - pos_rel - offset
  }
  
  # -------------------------------------------------------
  # MaxentScan per una sola variant
  # -------------------------------------------------------
  MaxentScanOne <- function(CHROMpos) {
    chr   <- paste0("chr", CHROMpos$CHROM)
    start <- CHROMpos$POS
    Ref   <- as.character(CHROMpos$REF)
    Alt   <- as.character(CHROMpos$ALT)
    query <- as.character(CHROMpos$genesList)
    gene  <- total_unique[match(query, total_unique[, 1]), 2]
    
    results <- list(
      acceptor_loss     = 0,  acceptor_loss_pos = NA,
      acceptor_gain     = 0,  acceptor_gain_pos = NA,
      donor_loss        = 0,  donor_loss_pos    = NA,
      donor_gain        = 0,  donor_gain_pos    = NA
    )
    
    pos1 <- as.numeric(start) - 26
    pos2 <- as.numeric(start) + 25 + nchar(Ref)
    
    if (!is.na(getStrand(gene))) {
      
      seq <- BSgenome::getSeq(BSgenome_assembly, chr,
                              start = as.integer(pos1),
                              end   = as.integer(pos2))
      
      pos1_alt <- as.integer(start) + nchar(Ref)
      pos2_alt <- pos1_alt + 26
      seq2     <- BSgenome::getSeq(BSgenome_assembly, chr,
                                   start = as.integer(pos1_alt),
                                   end   = as.integer(pos2_alt))
      
      seqALT <- Biostrings::DNAString(paste(seq[1:26], Alt, seq2, sep = ""))
      
      if (getStrand(gene) == "+") {
        r2a3 <- getBestScoreAndPos(Ref2Alt3(seq,    seqALT, Alt, Ref))
        r2a5 <- getBestScoreAndPos(Ref2Alt5(seq,    seqALT, Alt, Ref))
        a2r3 <- getBestScoreAndPos(Ref2Alt3(seqALT, seq,    Ref, Alt))
        a2r5 <- getBestScoreAndPos(Ref2Alt5(seqALT, seq,    Ref, Alt))
        
        results$acceptor_loss     <- -r2a3$score
        results$acceptor_loss_pos <- relToGenomic(r2a3$pos, start, "acceptor")
        results$acceptor_gain     <-  a2r3$score
        results$acceptor_gain_pos <- relToGenomic(a2r3$pos, start, "acceptor")
        results$donor_loss        <- -r2a5$score
        results$donor_loss_pos    <- relToGenomic(r2a5$pos, start, "donor")
        results$donor_gain        <-  a2r5$score
        results$donor_gain_pos    <- relToGenomic(a2r5$pos, start, "donor")
      }
      
      if (getStrand(gene) == "-") {
        seqComplement    <- Biostrings::reverseComplement(seq)
        seqALTComplement <- Biostrings::reverseComplement(seqALT)
        RefCompl <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(Ref)))
        AltCompl <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(Alt)))
        
        r2a3 <- getBestScoreAndPos(Ref2Alt3(seqComplement,    seqALTComplement, AltCompl, RefCompl))
        r2a5 <- getBestScoreAndPos(Ref2Alt5(seqComplement,    seqALTComplement, AltCompl, RefCompl))
        a2r3 <- getBestScoreAndPos(Ref2Alt3(seqALTComplement, seqComplement,    RefCompl, AltCompl))
        a2r5 <- getBestScoreAndPos(Ref2Alt5(seqALTComplement, seqComplement,    RefCompl, AltCompl))
        
        results$acceptor_loss     <- -r2a3$score
        results$acceptor_loss_pos <- relToGenomicMinus(r2a3$pos, start, "acceptor")
        results$acceptor_gain     <-  a2r3$score
        results$acceptor_gain_pos <- relToGenomicMinus(a2r3$pos, start, "acceptor")
        results$donor_loss        <- -r2a5$score
        results$donor_loss_pos    <- relToGenomicMinus(r2a5$pos, start, "donor")
        results$donor_gain        <-  a2r5$score
        results$donor_gain_pos    <- relToGenomicMinus(a2r5$pos, start, "donor")
      }
    }
    
    results$acceptor_loss <- round(results$acceptor_loss, 2)
    results$acceptor_gain <- round(results$acceptor_gain, 2)
    results$donor_loss    <- round(results$donor_loss,    2)
    results$donor_gain    <- round(results$donor_gain,    2)
    
    matrix(unlist(results), nrow = 1, ncol = 8)
  }
  
  # -------------------------------------------------------
  # Processa totes les variants del chunk
  # -------------------------------------------------------
  acumulat <- vector("list", nrow(matrixVariants_chunk))
  for (x in seq_len(nrow(matrixVariants_chunk))) {
    acumulat[[x]] <- cbind(matrixVariants_chunk[x, ], MaxentScanOne(matrixVariants_chunk[x, ]))
  }
  
  do.call(rbind, acumulat)
}


#' Evaluate splice site strength changes using MaxEntScan (parallel)
#'
#' @description
#' Main function for splice site evaluation in ClinPrior. Splits variants into
#' chunks and processes them in parallel using \code{BiocParallel}, calling
#' \code{.MaxentScanWorker()} in each worker. For each variant, computes
#' fractional MaxEntScan score changes for donor and acceptor splice sites
#' (loss and gain), taking strand into account. Typically called internally
#' by \code{priorBestVariant()} when \code{splicingMode = "full"}.
#'
#' @param matrixVariants A \code{data.frame} with columns \code{CHROM},
#'   \code{POS}, \code{REF}, \code{ALT}, and \code{genesList} (HGNC gene symbol).
#' @param assembly Character. Genome assembly version. Either \code{"assembly38"}
#'   (default, hg38) or \code{"assembly37"} (hg19).
#' @param Threshold Numeric. Minimum fractional score change to consider a
#'   splice site as affected. Default \code{0.2}.
#' @param chunk_size Integer. Number of variants processed per parallel chunk.
#'   Default \code{250}.
#' @param workers Integer. Number of parallel workers. Defaults to the number
#'   of available cores minus one. Use \code{1} for sequential (debug) mode.
#' @param BPPARAM A \code{BiocParallelParam} object. If \code{NULL} (default),
#'   a \code{SnowParam} (or \code{SerialParam} if \code{workers = 1}) is created
#'   and destroyed locally. Pass an external \code{BatchtoolsParam} for HPC/SLURM.
#'
#' @return A \code{data.frame} with the original \code{matrixVariants} columns
#'   plus 8 additional columns: \code{acceptor_loss}, \code{acceptor_loss_pos},
#'   \code{acceptor_gain}, \code{acceptor_gain_pos}, \code{donor_loss},
#'   \code{donor_loss_pos}, \code{donor_gain}, \code{donor_gain_pos}.
#'   Score columns contain the fractional change (0 = no effect). Position
#'   columns contain the genomic coordinate of the affected splice site
#'   (\code{NA} if no effect detected).
#'
#' @export
#' @importFrom BiocParallel bplapply SnowParam SerialParam bpisup bpstart bpstop bpnworkers
#'
#' @examples
#' \dontrun{
#' # Called automatically by priorBestVariant() with splicingMode = "full".
#' # Can also be run standalone:
#' mx_input <- data.frame(CHROM = "1", POS = 123456, REF = "A", ALT = "T",
#'                        genesList = "BRCA2")
#' output <- MaxentScanClinPrior(mx_input, assembly = "assembly38", workers = 4)
#' }
MaxentScanClinPrior <- function(matrixVariants,
                                assembly   = "assembly38",
                                Threshold  = 0.2,
                                chunk_size = 250,
                                workers    = max(1, parallel::detectCores() - 1),
                                BPPARAM    = NULL) {
  
  # --Gestio del cicle de vida del BPPARAM --------------------------------
  external_bpparam <- !is.null(BPPARAM)
  
  if (!external_bpparam) {
    if (workers == 1) {
      BPPARAM <- BiocParallel::SerialParam()
    } else {
      BPPARAM <- BiocParallel::SnowParam(workers = workers, progressbar = TRUE)
    }
  }
  
  on.exit({
    if (!external_bpparam) {
      tryCatch(
        { if (BiocParallel::bpisup(BPPARAM)) BiocParallel::bpstop(BPPARAM) },
        error = function(e) NULL
      )
    }
    gc()
  }, add = TRUE)
  
  if (!BiocParallel::bpisup(BPPARAM)) BiocParallel::bpstart(BPPARAM)
  
  maxentpyPath <- system.file("extdata", package = "ClinPrior")
  
  n      <- nrow(matrixVariants)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  
  message(sprintf(
    "Processant %d variants en %d chunks (chunk_size = %d) | backend: %s | workers: %d",
    n, length(chunks), chunk_size, class(BPPARAM), BiocParallel::bpnworkers(BPPARAM)
  ))
  
  results <- BiocParallel::bplapply(
    chunks,
    FUN            = .MaxentScanWorker,
    matrixVariants = matrixVariants,
    assembly       = assembly,
    Threshold      = Threshold,
    maxentpyPath   = maxentpyPath,
    BPPARAM        = BPPARAM
  )
  
  output <- do.call(rbind, results)
  
  col_start <- ncol(matrixVariants) + 1
  colnames(output)[col_start:(col_start + 7)] <- c(
    "acceptor_loss", "acceptor_loss_pos",
    "acceptor_gain", "acceptor_gain_pos",
    "donor_loss",    "donor_loss_pos",
    "donor_gain",    "donor_gain_pos"
  )
  
  return(output)
}