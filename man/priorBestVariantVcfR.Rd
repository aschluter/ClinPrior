% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priorBestVariantsVcfR.r
\name{priorBestVariantVcfR}
\alias{priorBestVariantVcfR}
\title{Title}
\usage{
priorBestVariantVcfR(
  variants,
  sampleName,
  filter = "",
  isCodingVar = TRUE,
  frequenceAR = 0.01,
  frequenceAD = 5e-05,
  GlobalPhenotypicScore,
  assembly = "assembly37",
  processors = 1
)
}
\arguments{
\item{variants}{object of class vcfR-class with the patient variants.}

\item{sampleName}{character name with the patient code as is written in the VCF file.}

\item{filter}{character name with the desired filter to select the variants. It should be as is written in the FILTER column in the vcf file.}

\item{isCodingVar}{logical indicating whether to consider only variants in coding regions. Default TRUE.}

\item{frequenceAR}{numerical value indicating the MAF(minor allele frequence) to filter out the variants in autosomal recessive and X-linked inheritance. Default 0.01}

\item{frequenceAD}{numerical value indicating the MAF(minor allele frequence) to filter out the variants in autosomal dominant recessive inheritance. Default 0.00005}

\item{GlobalPhenotypicScore}{matrix with the phenotypic metrics obtained from the MatrixPropagation function.}

\item{assembly}{Genome assembly used. Default assembly human GRCh37.}

\item{processors}{Default=1.}
}
\value{
data frame with the patient's variants and their associated information classified from most to least likely to be the cause of the patient's phenotype.
}
\description{
Title
}
\examples{
library(vcfR)
patientHPOsFile <- paste(system.file("extdata/example", package = "ClinPrior"),"HPOpatient.txt",sep="/")
HPOpatient <- unique(read.csv(patientHPOsFile, header = FALSE, sep = "\t")[, 1])
Y<-proteinScore(HPOpatient)
ClinPriorGeneScore<-MatrixPropagation(Y,alpha=0.2)

vcfFile = paste(system.file("extdata/example", package = "ClinPrior"),"HG001_GRCh37_1_22_v4.2.1_benchmark.vep01.KCNQ2Met546Thr.vcf.gz",sep="/")
variants <- read.vcfR(vcfFile)
variantsFiltered <- readVCF(sampleName = "HG001",variants=variants)
result = priorBestVariantVcfR(variants = variantsFiltered, sampleName = "HG001",GlobalPhenotypicScore = ClinPriorGeneScore)
}
