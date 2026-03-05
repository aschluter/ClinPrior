# Suppress R CMD check notes for dataset variables used across functions
utils::globalVariables(c(
  "totalEnsembl", "geneConstraints", "HPO2genes", "HPOadj",
  "HPOdistance", "HPOqueryGene", "HPOorigGroups", "treureHPO",
  "total_unique", "total_unique_Conn_Physical", "total_unique_Conn_Func",
  "normPhysical", "normFunc", "gene2strand", "g1",
  "hg38.omim2gene.approvedGeneSymbol", "hg38.omimPhenotype.phenotypeId",
  "hg38.omimGeneMap2.phenotypes", "approvedGeneSymbol", "phenotype",
  "phenotypeId", "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38",
  "chr", "start", "end", "z", "idx", "vp", "vc", "HPOorigGroups",
  "py", "import_from_path", "py_run_string"
))
utils::globalVariables(c(
  # ... les que ja tenies ...
  "."
))
