# ClinPrior data update
#
# Rebuilds data/ClinPrior.RData with the latest versions of:
#   - HPO ontology (hp.obo)                    -> g1, HPOdistance
#   - HPO gene annotations (genes_to_phenotype) -> HPO2genes, HPOadj, HPOqueryGene
#   - gnomAD constraint metrics                -> geneConstraints
#   - Ensembl BioMart (archive 116)            -> totalEnsembl, gene2strand
#   - optional: BioGRID multi-validated physical interactions and BioPlex
#     networks (HEK293T 3.0, HCT116), added to the physical network
#     (normPhysical, total_unique_Conn_Physical)
# The functional network (normFunc, total_unique_Conn_Func), treureHPO and,
# without --biogrid/--bioplex, the physical network are kept from the
# current RData.
# Genes annotated in HPO that are not in total_unique are added to it.
#
# Usage (from the package root):
#   Rscript data-raw/update_clinprior.R [options]
# Options:
#   --rdata=PATH       RData to update (default data/ClinPrior.RData)
#   --work=DIR         download folder (default a temporary folder)
#   --report=PATH      Markdown report (default <work>/update_report.md)
#   --gnomad-url=URL   gnomAD constraint table
#   --skip-gnomad      keep geneConstraints from the current RData
#   --skip-ensembl     keep totalEnsembl and gene2strand from the current RData
#                      (also kept, with a note, if Ensembl does not answer)
#   --ensembl-host=URL Ensembl BioMart host (default https://e116.ensembl.org)
#   --biogrid          add BioGRID multi-validated (MV) physical interactions
#   --biogrid-systems=binary|ppi
#                      binary (default): direct protein-protein systems;
#                      ppi: all physical systems except co-localization and RNA
#   --biogrid-url=URL  BioGRID MV physical tab3 zip (default: latest release)
#   --bioplex[=293T,HCT116]
#                      add BioPlex networks (default: both cell lines)
#   --bioplex-min-pint=X
#                      minimum BioPlex interaction probability (default 0:
#                      keep the published networks as they are)
#   --dry-run          write the report but not the RData
# Exit status: 0 = RData updated, 3 = no changes, other = error.

suppressPackageStartupMessages({
  library(Matrix)
  library(igraph)
})

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "(=|$)"), args, value = TRUE)
  if(length(hit) == 0) return(default)
  if(!grepl("=", hit[1])) return(TRUE)
  sub(paste0("^--", name, "="), "", hit[1])
}

rdata      <- opt("rdata", "data/ClinPrior.RData")
work       <- opt("work", file.path(tempdir(), "clinprior_update"))
report     <- opt("report", file.path(work, "update_report.md"))
gnomad_url <- opt("gnomad-url",
                  "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv")
skip_gnomad  <- isTRUE(opt("skip-gnomad", FALSE))
skip_ensembl <- isTRUE(opt("skip-ensembl", FALSE))
dry_run      <- isTRUE(opt("dry-run", FALSE))
use_biogrid  <- isTRUE(opt("biogrid", FALSE))
bioplex_opt  <- opt("bioplex", FALSE)
use_bioplex  <- !isFALSE(bioplex_opt)
bioplex_cl   <- if(isTRUE(bioplex_opt)) c("293T", "HCT116") else strsplit(as.character(bioplex_opt), ",")[[1]]
bioplex_min  <- as.numeric(opt("bioplex-min-pint", 0))
biogrid_sys  <- opt("biogrid-systems", "binary")
ensembl_host <- opt("ensembl-host", "https://e116.ensembl.org")
biogrid_url  <- opt("biogrid-url",
                    "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip")

hpo_url <- "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/"

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
for(f in list.files(file.path(here, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

dir.create(work, recursive = TRUE, showWarnings = FALSE)
download <- function(url, file) {
  dest <- file.path(work, file)
  if(!file.exists(dest)) {
    message("Downloading ", url)
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

# -- Current data -------------------------------------------------------------
old <- new.env()
load(rdata, envir = old)
new <- new.env()
for(n in ls(old)) assign(n, get(n, old), envir = new)
log <- character(0)
note <- function(...) { txt <- paste0(...); message(txt); log <<- c(log, txt) }

# -- HPO ----------------------------------------------------------------------
obo <- download(paste0(hpo_url, "hp.obo"), "hp.obo")
g2p <- download(paste0(hpo_url, "genes_to_phenotype.txt"), "genes_to_phenotype.txt")
hpo_version <- sub("^data-version: ", "",
                   grep("^data-version:", readLines(obo, n = 30), value = TRUE)[1])

new$g1        <- buildHPOgraph(obo)
new$HPO2genes <- updateHPO2genes(g2p, old$HPO2genes)
symbols       <- attr(new$HPO2genes, "symbols")
attr(new$HPO2genes, "symbols") <- NULL

# Genes with HPO annotation must be in total_unique (proteinScore() needs it)
missing <- setdiff(unique(new$HPO2genes$V1), new$total_unique[, 2])
if(length(missing) > 0) {
  sym <- symbols$gene_symbol[match(missing, symbols$ncbi_gene_id)]
  sym[is.na(sym)] <- as.character(missing[is.na(sym)])
  add <- data.frame(missing, missing, stringsAsFactors = FALSE)
  names(add) <- names(new$total_unique)
  add[, 1] <- sym
  new$total_unique <- rbind(new$total_unique, add)
  rownames(new$total_unique) <- NULL
  note("Genes added to total_unique (annotated in HPO, not in the gene list): ",
       paste0(sym, " (", missing, ")", collapse = ", "))
}

adj <- buildHPOadj(new$HPO2genes, new$g1, new$treureHPO)
for(n in names(adj)) assign(n, adj[[n]], envir = new)

# -- Physical network additions (optional): BioGRID MV and BioPlex -----------
biogrid_release <- NA
bioplex_release <- NA
extra_pairs <- list()
if(use_biogrid) {
  zip <- download(biogrid_url, "BIOGRID-MV-Physical.tab3.zip")
  extra_pairs$BioGRID <- readBioGRID(zip, biogrid_sys)
  biogrid_release <- attr(extra_pairs$BioGRID, "release")
  note("BioGRID (", biogrid_release, ", systems: ", biogrid_sys, "): ",
       nrow(extra_pairs$BioGRID), " human interactions")
}
if(use_bioplex) {
  bad <- setdiff(bioplex_cl, names(bioplexUrls))
  if(length(bad)) stop("Unknown BioPlex cell line(s): ", paste(bad, collapse = ", "))
  files <- vapply(bioplex_cl, function(cl) download(bioplexUrls[[cl]], basename(bioplexUrls[[cl]])), "")
  extra_pairs$BioPlex <- readBioPlex(files, bioplex_min)
  bioplex_release <- attr(extra_pairs$BioPlex, "release")
  note("BioPlex (", bioplex_release, if(bioplex_min > 0) paste0(", pInt >= ", bioplex_min) else "",
       "): ", nrow(extra_pairs$BioPlex), " interactions")
}
if(length(extra_pairs) > 0) {
  pairs <- unique(do.call(rbind, lapply(extra_pairs, function(p) p[, c("a", "b")])))
  if(length(extra_pairs) == 2) {
    k <- function(p) paste(p$a, p$b)
    note("Interactions in both BioGRID and BioPlex: ",
         length(intersect(k(extra_pairs$BioGRID), k(extra_pairs$BioPlex))))
  }
  net <- addEdgesToNetwork(new$normPhysical, new$total_unique_Conn_Physical, pairs)
  new$normPhysical <- net$W
  new$total_unique_Conn_Physical <- as.integer(net$genes)
  note("Physical network: ", net$added_edges, " new edges, ", net$present_edges,
       " already present; ", length(net$new_genes), " genes added",
       if(length(net$capped)) paste0("; no new edges for genes with more than 1000 neighbours: ",
                                     paste(net$capped, collapse = ", ")) else "")
  # Network genes must be in total_unique
  extra <- setdiff(net$genes, new$total_unique[, 2])
  if(length(extra) > 0) {
    sym <- do.call(rbind, lapply(extra_pairs, attr, "symbols"))
    s2  <- sym$symbol[match(extra, sym$entrez)]
    s2[is.na(s2) | s2 == ""] <- as.character(extra[is.na(s2) | s2 == ""])
    add <- data.frame(s2, as.integer(extra), stringsAsFactors = FALSE)
    names(add) <- names(new$total_unique)
    new$total_unique <- rbind(new$total_unique, add)
    rownames(new$total_unique) <- NULL
    note("Genes added to total_unique from the physical network sources: ", length(extra))
  }
}

# -- gnomAD -------------------------------------------------------------------
if(!skip_gnomad) {
  tsv <- download(gnomad_url, basename(gnomad_url))
  new$geneConstraints <- buildGeneConstraints(tsv)
} else note("geneConstraints kept from the current RData (--skip-gnomad)")

# -- Ensembl ------------------------------------------------------------------
ensembl_version <- NA
if(!skip_ensembl) {
  bm <- tryCatch(ensemblBioMart(ensembl_host), error = function(e) e)
  if(inherits(bm, "error")) {
    note("Ensembl BioMart did not answer (", conditionMessage(bm),
         "): totalEnsembl and gene2strand kept from the current RData")
  } else {
    rs <- refreshSymbols(new$total_unique, bm)
    new$total_unique <- rs$total_unique
    if(nrow(rs$changes) > 0) {
      write.table(rs$changes, file.path(work, "symbol_changes.tsv"), sep = "\t",
                  quote = FALSE, row.names = FALSE)
      ex <- head(paste(rs$changes$old, "->", rs$changes$new), 30)
      note("Symbols updated in total_unique: ", nrow(rs$changes), " (",
           paste(ex, collapse = ", "), if(nrow(rs$changes) > 30) ", ...", ")")
    }
    new$totalEnsembl <- buildTotalEnsembl(new$total_unique, bm)
    new$gene2strand  <- buildGene2strand(new$total_unique, bm)
    ensembl_version  <- ensembl_host
  }
} else note("totalEnsembl and gene2strand kept from the current RData (--skip-ensembl)")

# -- Checks -------------------------------------------------------------------
check <- function(ok, msg) if(!isTRUE(ok)) stop("Check failed: ", msg, call. = FALSE)
check(all(unique(new$HPO2genes$V1) %in% new$total_unique[, 2]), "HPO2genes genes not in total_unique")
check(identical(colnames(new$HPOadj), V(new$g1)$name), "HPOadj columns != g1 vertices")
check(nrow(new$HPOqueryGene) == nrow(new$HPOadj), "HPOqueryGene rows != HPOadj rows")
check(all(c(new$total_unique_Conn_Physical, new$total_unique_Conn_Func) %in% new$total_unique[, 2]),
      "network genes not in total_unique")
for(n in c("HPO2genes", "geneConstraints", "totalEnsembl", "gene2strand")) {
  check(NROW(get(n, new)) >= 0.9 * NROW(get(n, old)),
        paste0(n, " has less than 90% of the previous rows (", NROW(get(n, new)), " vs ", NROW(get(n, old)), ")"))
}

# -- Compare and save ---------------------------------------------------------
objs    <- c("g1", "treureHPO", "HPOadj", "HPOqueryGene", "gene2strand", "HPOdistance",
             "geneConstraints", "HPO2genes", "normPhysical", "normFunc", "total_unique",
             "totalEnsembl", "total_unique_Conn_Physical", "total_unique_Conn_Func")
size    <- function(o) paste(if(is.null(dim(o))) length(o) else dim(o), collapse = " x ")
same    <- function(n) {
  a <- get(n, old); b <- get(n, new)
  if(inherits(a, "igraph")) return(identical(V(a)$name, V(b)$name) && identical(as_edgelist(a), as_edgelist(b)))
  identical(a, b)
}
changed <- objs[!vapply(objs, same, TRUE)]

rows <- vapply(objs, function(n) sprintf("| `%s` | %s | %s | %s |", n, size(get(n, old)),
                                         size(get(n, new)), if(n %in% changed) "**updated**" else "="), "")
gene_add <- length(setdiff(unique(new$HPO2genes$V1), unique(old$HPO2genes$V1)))
pair_add <- nrow(new$HPO2genes) - nrow(old$HPO2genes)
md <- c(
  "# ClinPrior data update", "",
  paste0("Date: ", format(Sys.Date())), "",
  "## Sources", "",
  paste0("- HPO: `", hpo_version, "` (hp.obo, genes_to_phenotype.txt)"),
  paste0("- gnomAD constraint: ", if(skip_gnomad) "not updated" else paste0("`", basename(gnomad_url), "`")),
  paste0("- Ensembl BioMart: ", if(is.na(ensembl_version)) "not updated" else ensembl_version),
  paste0("- BioGRID MV physical: ", if(use_biogrid) paste0("`", biogrid_release, "` (", biogrid_sys, ")") else "not used"),
  paste0("- BioPlex: ", if(use_bioplex) bioplex_release else "not used"),
  "",
  "## Objects", "",
  "| Object | Before | After | |", "|---|---|---|---|", rows, "",
  paste0("HPO2genes: ", pair_add, " new gene-HPO pairs, ", gene_add, " new genes."),
  if(length(log)) c("", "## Notes", "", paste0("- ", log)) else character(0)
)
dir.create(dirname(report), recursive = TRUE, showWarnings = FALSE)
writeLines(md, report)
cat(md, sep = "\n")

if(length(changed) == 0) {
  message("No changes: ", rdata, " not written")
  quit(status = 3)
}
if(dry_run) {
  message("--dry-run: ", rdata, " not written")
  quit(status = 0)
}
save(list = objs, envir = new, file = rdata, compress = "xz")
message("Written ", rdata, " (updated: ", paste(changed, collapse = ", "), ")")