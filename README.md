

# ClinPrior

ClinPrior is an interactome-driven prioritization method that predicts the patient's disease gene based on the description of the phenotype in HPO terms. This prioritization is divided into two steps: (1) computation of a phenotypic metric by comparing the patient's phenotype with HPO-Gene associations from existing human disease databases (prior knowledge) and (2) iterative propagation of this phenotypic score within of a multilayer network with physical and functional interactions.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("aschluter/ClinPrior")
```


## Gene prioritization

Gene prioritization according HPOs provided by user. 

``` r
Y<-proteinScore(<HPOpatient>)
```


## Score propagation

Propagation of the phenotypic score in the physical and functional multilayer network. 

``` r
PriorPhenoScore<-MatrixPropagation(<Y>,alpha=0.2)
```


