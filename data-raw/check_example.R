# Regression check after a data update: the package example (KCNQ2
# p.Met546Thr, homozygous) must still be ranked first in both assemblies.
#
# Usage (from the package root): Rscript data-raw/check_example.R [report.md]
# Appends the result to the report if given. Exit status 1 if the check fails.

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))

report <- commandArgs(trailingOnly = TRUE)[1]
ex  <- system.file("extdata/example", package = "ClinPrior")
hpo <- readLines(file.path(ex, "HPOpatient.txt"))
hpo <- hpo[nzchar(hpo)]
G   <- suppressMessages(MatrixPropagation(proteinScore(hpo)$scores, alpha = 0.2))

res <- lapply(c("37", "38"), function(a) {
  vcf <- file.path(ex, paste0("HG001_GRCh", a, "_1_22_v4.2.1_benchmark.vep01.KCNQ2Met546Thr.vcf.gz"))
  v <- suppressMessages(readVCF(sampleName = "HG001", vcfFile = vcf, assembly = paste0("assembly", a)))
  r <- suppressMessages(priorBestVariant(v, sampleName = "HG001", GlobalPhenotypicScore = G,
                                         assembly = paste0("assembly", a), splicingMode = "vep",
                                         discard_zero_score = TRUE))
  k <- which(r$genesList == "KCNQ2")[1]
  data.frame(assembly = paste0("GRCh", a), variants = nrow(r),
             KCNQ2_rank = r$ClinPriorPosition[k], KCNQ2_score = r$ClinPriorScore[k],
             first = r$genesList[1])
})
res <- do.call(rbind, res)
print(res, row.names = FALSE)

ok <- all(res$KCNQ2_rank %in% 1)
md <- c("", "## Example check", "",
        "| Assembly | Variants | KCNQ2 rank | KCNQ2 score |", "|---|---|---|---|",
        sprintf("| %s | %d | %s | %s |", res$assembly, res$variants, res$KCNQ2_rank, res$KCNQ2_score),
        "", if(ok) "KCNQ2 is ranked first in both assemblies." else
          "**KCNQ2 is not ranked first: review the update before merging.**")
if(!is.na(report)) cat(md, file = report, sep = "\n", append = TRUE)
if(!ok) quit(status = 1)