# Gene constraint metrics from the gnomAD constraint table
# (gnomad.v4.1.constraint_metrics.tsv), canonical transcripts only.
# gnomAD v4.1 has Ensembl and RefSeq transcripts, so "gene" is an Ensembl
# gene ID or an Entrez ID; priorBestVariant() looks up both.
buildGeneConstraints <- function(gnomad_tsv) {
  g <- utils::read.delim(gnomad_tsv, sep = "\t", header = TRUE, quote = "",
                         stringsAsFactors = FALSE)
  g <- g[tolower(as.character(g$canonical)) == "true", ]
  out <- data.frame(gene  = as.character(g$gene_id),
                    lof_z = g$lof.z_score, mis_z = g$mis.z_score,
                    pLI   = g$lof_hc_lc.pLI, pRec = g$lof_hc_lc.pRec,
                    pNull = g$lof_hc_lc.pNull, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}