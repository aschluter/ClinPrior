# Update HPO2genes with the HPO gene annotations (genes_to_phenotype.txt).
#
# As in the original update, new gene-HPO pairs are added to the previous
# table (pairs are never removed) and rows are ordered by the Entrez ID as
# text, which defines the row order of HPOadj.
updateHPO2genes <- function(genes_to_phenotype, old = NULL) {
  g2p <- utils::read.delim(genes_to_phenotype, sep = "\t", header = TRUE,
                           quote = "", stringsAsFactors = FALSE)
  new <- unique(data.frame(V1 = as.integer(g2p$ncbi_gene_id), V2 = g2p$hpo_id,
                           stringsAsFactors = FALSE))
  new <- new[!is.na(new$V1) & grepl("^HP:[0-9]{7}$", new$V2), ]
  new <- new[order(as.character(new$V1), method = "radix"), ]
  out <- if(is.null(old)) new else
    unique(rbind(data.frame(V1 = as.integer(old[, 1]), V2 = old[, 2], stringsAsFactors = FALSE), new))
  out <- out[order(as.character(out$V1), method = "radix"), ]
  rownames(out) <- NULL
  attr(out, "symbols") <- unique(g2p[, c("ncbi_gene_id", "gene_symbol")])
  out
}