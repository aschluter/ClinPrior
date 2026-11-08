# Gene identifiers and strand from Ensembl BioMart, for the Entrez IDs of
# total_unique.
#
# Ensembl 116 (June 2026) is the last release of the legacy Ensembl site,
# where BioMart runs; www.ensembl.org is now the new site, without BioMart,
# and the USEast/Asia mirrors are retired. The query is sent directly to the
# archive host (HTTP GET), without biomaRt.

# All human genes: hgnc_symbol, entrezgene_id, ensembl_gene_id, strand
ensemblBioMart <- function(host = "https://e116.ensembl.org", tries = 3) {
  q <- paste0('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>',
              '<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" ',
              'completionStamp="1" datasetConfigVersion="0.6">',
              '<Dataset name="hsapiens_gene_ensembl" interface="default">',
              '<Attribute name="hgnc_symbol"/><Attribute name="entrezgene_id"/>',
              '<Attribute name="ensembl_gene_id"/><Attribute name="strand"/>',
              '</Dataset></Query>')
  url <- paste0(host, "/biomart/martservice?query=", utils::URLencode(q, reserved = TRUE))
  for(k in seq_len(tries)) {
    r <- tryCatch(curl::curl_fetch_memory(url, handle = curl::new_handle(timeout = 900)),
                  error = function(e) e)
    if(!inherits(r, "error") && r$status_code == 200) {
      lines <- strsplit(rawToChar(r$content), "\n", fixed = TRUE)[[1]]
      # A complete answer starts with the header and ends with "[success]"
      if(length(lines) > 2 && startsWith(lines[1], "HGNC symbol") &&
         utils::tail(lines, 1) == "[success]") {
        x <- utils::read.delim(text = lines[-c(1, length(lines))], header = FALSE,
                               quote = "", fill = TRUE, colClasses = "character",
                               col.names = c("hgnc_symbol", "entrezgene_id",
                                             "ensembl_gene_id", "strand"))
        attr(x, "host") <- host
        return(x)
      }
      msg <- paste("unexpected answer:", substr(paste(lines[1:min(2, length(lines))], collapse = " "), 1, 200))
    } else {
      msg <- if(inherits(r, "error")) conditionMessage(r) else paste("HTTP", r$status_code)
    }
    message("Ensembl BioMart (", host, ") not available: ", msg, if(k < tries) "; retry in 60 s")
    if(k < tries) Sys.sleep(60)
  }
  stop("Could not get the gene table from Ensembl BioMart (", host, ")")
}

# totalEnsembl: hgnc_symbol, entrezgene_id, ensembl_gene_id
buildTotalEnsembl <- function(total_unique, bm = ensemblBioMart()) {
  x <- bm[bm$entrezgene_id != "" & bm$entrezgene_id %in% as.character(total_unique[, 2]), ]
  out <- unique(data.frame(hgnc_symbol = x$hgnc_symbol,
                           entrezgene_id = as.integer(x$entrezgene_id),
                           ensembl_gene_id = x$ensembl_gene_id, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}

# gene2strand: Entrez ID (column named "symbol", as used by the package) and +/-
buildGene2strand <- function(total_unique, bm = ensemblBioMart()) {
  x <- bm[bm$entrezgene_id != "" & bm$entrezgene_id %in% as.character(total_unique[, 2]), ]
  x <- unique(data.frame(symbol = as.integer(x$entrezgene_id), hgnc = x$hgnc_symbol,
                         strand = ifelse(x$strand == "1", "+", "-"), stringsAsFactors = FALSE))
  # An Entrez ID linked to several Ensembl genes: keep the gene with the same
  # symbol as total_unique (one row per Entrez ID)
  same <- x$hgnc == total_unique[match(x$symbol, total_unique[, 2]), 1]
  same[is.na(same)] <- FALSE
  x   <- x[order(x$symbol, !same), ]
  out <- x[!duplicated(x$symbol), c("symbol", "strand")]
  out <- out[order(as.character(out$symbol), method = "radix"), ]
  rownames(out) <- NULL
  out
}

# Refresh the symbols of total_unique (column 1) with the HGNC symbols of
# Ensembl. Rows and Entrez IDs (column 2) are not changed. A symbol is only
# replaced when Ensembl links the Entrez ID to a single HGNC symbol and that
# symbol is not used by another gene; otherwise the current one is kept.
refreshSymbols <- function(total_unique, bm = ensemblBioMart()) {
  x <- unique(bm[bm$entrezgene_id != "" & bm$hgnc_symbol != "", c("entrezgene_id", "hgnc_symbol")])
  n <- table(x$entrezgene_id)
  x <- x[x$entrezgene_id %in% names(n)[n == 1], ]
  cur <- as.character(total_unique[, 1])
  sym <- x$hgnc_symbol[match(as.character(total_unique[, 2]), x$entrezgene_id)]
  upd <- !is.na(sym) & sym != cur
  repeat {                       # never give one symbol to two genes
    final <- ifelse(upd, sym, cur)
    clash <- upd & final %in% final[duplicated(final)]
    if(!any(clash)) break
    upd[clash] <- FALSE
  }
  out <- total_unique
  out[upd, 1] <- sym[upd]
  list(total_unique = out,
       changes = data.frame(geneID = total_unique[upd, 2], old = cur[upd], new = sym[upd]))
}