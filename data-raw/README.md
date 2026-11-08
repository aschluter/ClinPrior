# ClinPrior data update

Scripts to rebuild `data/ClinPrior.RData` with new releases of the source
databases. They are not part of the installed package (`data-raw` is in
`.Rbuildignore`).

## What is updated

| Object | Source | Script |
|---|---|---|
| `g1` | HPO ontology, [`hp.obo`](https://github.com/obophenotype/human-phenotype-ontology/releases) | `R/hpo_graph.R` |
| `HPO2genes` | HPO gene annotations, `genes_to_phenotype.txt` | `R/hpo2genes.R` |
| `HPOadj`, `HPOdistance`, `HPOqueryGene` | `HPO2genes` + `g1` | `R/hpo_adj.R` |
| `geneConstraints` | [gnomAD v4.1 constraint metrics](https://gnomad.broadinstitute.org/downloads#v4-constraint) | `R/gene_constraints.R` |
| `totalEnsembl`, `gene2strand` | Ensembl BioMart, archive of Ensembl 116 | `R/ensembl.R` |
| `normPhysical`, `total_unique_Conn_Physical` (optional) | [BioGRID MV physical](https://wiki.thebiogrid.org/doku.php/biogrid_mv), [BioPlex](https://bioplex.hms.harvard.edu/interactions.php) | `R/biogrid.R`, `R/bioplex.R` |

Kept from the current RData: `normFunc` and `total_unique_Conn_Func`
(functional network), `treureHPO` and, unless `--biogrid` or `--bioplex`
is used, `normPhysical` and `total_unique_Conn_Physical` (physical
network). `total_unique` is kept and extended with the genes annotated in
HPO that are not yet in it (`proteinScore()` needs every gene of
`HPO2genes` in `total_unique`). Its rows and Entrez IDs never change, but
its symbols are refreshed with the HGNC symbols of Ensembl (when Ensembl
links the Entrez ID to a single symbol not used by another gene), so that
they match `totalEnsembl`: `MaxentScanClinPrior()` looks up genes by symbol.

Notes on the rules:

- `HPO2genes` accumulates: new gene-HPO pairs are added to the previous
  table and no pair is removed.
- `gene2strand` has one row per Entrez ID (`MaxentScanClinPrior()` takes
  the first match). When an Entrez ID is linked to several Ensembl genes,
  the gene with the same symbol as `total_unique` is used.
- Row order matters: `HPOadj` and `HPOqueryGene` rows follow the genes of
  `HPO2genes` (after removing `treureHPO`) in table order, and
  `HPOadj`/`HPOdistance` columns follow the vertices of `g1`.
  `proteinScore()` relies on both.
- `HPOdistance` is 1 when two HPO terms are the same or directly linked in
  `g1` (undirected distance < 2). It is built from the sparse adjacency
  matrix, which gives the same matrix as `igraph::distances()` without the
  dense HPO x HPO matrix (about 4 GB).
- `geneConstraints$gene` is an Ensembl gene ID or an Entrez ID (gnomAD v4.1
  has Ensembl and RefSeq transcripts); `priorBestVariant()` looks up both.
- Ensembl: Ensembl 116 (June 2026) is the last release of the legacy site,
  the only one with BioMart (www.ensembl.org is now the new site, without
  BioMart, and the USEast/Asia mirrors are retired). The table is queried
  directly from `https://e116.ensembl.org` (HTTP GET, without `biomaRt`);
  `--ensembl-host` changes the archive. If Ensembl does not answer, the
  previous `totalEnsembl` and `gene2strand` are kept and the report says so.
  For releases after 116, these tables will have to be built from the
  Ensembl FTP files.

## Gene networks

The networks are described in the ClinPrior paper (Genome Medicine 2023):
physical layer from BioPlex 2.0, HuRI, Lit-BM-13, HI-I-05 and HI-II-14;
functional layer from HumanNet-CF v2, KEGG, Recon3D and Signor 2.0. Only
the normalised matrices are stored: `norm = D^-1 (A + I)`, where `A` has
integer, symmetric edge weights (number of times the edge appears in the
sources, usually 2), `I` is a self loop of weight 1 and `D` the row sums.
`A` can therefore be recovered exactly (`R/biogrid.R`).

### Optional: BioGRID multi-validated interactions (`--biogrid`)

Adds to the physical layer the human interactions of the
[BioGRID multi-validated (MV) physical dataset](https://wiki.thebiogrid.org/doku.php/biogrid_mv)
(each interaction supported by at least two experimental systems or two
publications). By default (`--biogrid-systems=binary`) only direct
(binary) systems are used: two-hybrid, PCA, reconstituted complex,
co-crystal structure, FRET, far western, protein-peptide and biochemical
activity, as the binary sources of the original network (HuRI, HI,
Lit-BM-13). `ppi` adds co-complex systems from the literature (affinity
capture, co-fractionation, proximity labelling...), except
co-localization and RNA interactions; co-complex interactions are taken
from BioPlex instead (systematic screen, less biased towards well-studied
genes).

### Optional: BioPlex networks (`--bioplex`)

Adds the BioPlex affinity purification-mass spectrometry networks as
published by BioPlex: BioPlex 3.0 in HEK293T cells (which contains BioPlex
2.0, already in the original network) and BioPlex HCT116. `--bioplex=293T`
or `--bioplex=HCT116` selects one cell line; `--bioplex-min-pint=X` keeps
only interactions with interaction probability `pInt >= X` (by default the
published networks are used as they are).

BioGRID and BioPlex can be used together (recommended:
`--biogrid --bioplex`): their interactions are merged first, so an
interaction found in both sources is added once.

Rules for the added interactions: edges already in the network keep their
weight; new edges get weight 2; new genes are added to the network and to
`total_unique`; genes of the original network are never removed. Genes
that would end up with more than 1000 neighbours (the limit of the
original network) get no new edges and keep their original ones; new
genes in that situation are not added. The report lists them.

### Evaluation

Adding interactions changes the propagation for every patient. Evaluate
it on cases with a known diagnosis before using it in the monthly update.
It is off by default; in GitHub, it can be enabled when running the
workflow manually.

## Run locally

From the package root:

```bash
Rscript data-raw/update_clinprior.R            # downloads sources, updates data/ClinPrior.RData
Rscript data-raw/check_example.R               # KCNQ2 must be ranked first (GRCh37 and GRCh38)
```

Options of `update_clinprior.R`: `--work=DIR` (keep the downloads),
`--report=FILE`, `--skip-gnomad`, `--skip-ensembl`, `--ensembl-host=URL`,
`--dry-run`, `--gnomad-url=URL`, `--biogrid`, `--biogrid-systems=binary|ppi`,
`--biogrid-url=URL`, `--bioplex[=293T,HCT116]`, `--bioplex-min-pint=X`.
Exit status 3 means that nothing changed.

Requirements: R packages `Matrix`, `igraph`, `curl`, `pkgload`, and
`bcftools` for the example check.

## Automatic update

`.github/workflows/update-clinprior-data.yml` runs on day 1 of each month
(and manually from the Actions tab): it checks out `dev`, runs both scripts
and, if the data changed, opens a pull request to `dev` with the report.
Review the report and the example check before merging.