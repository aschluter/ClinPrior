

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

| HP:0004481 | Progressive macrocephaly |
| HP:0002376 | Developmental regression |
| HP:0001257 | Spasticity |
| HP:0001250 | Seizure |
| HP:0000238 | Hydrocephalus |
| HP:0002922 | Increased CSF protein |
| HP:0000365 | Hearing impairment |

``` r

HPOpatient = c("HP:0004481","HP:0002376","HP:0001257","HP:0001250","HP:0000238","HP:0002922","HP:0000365"")

Y<-proteinScore(HPOpatient)
```


## Score propagation

Propagation of the phenotypic score in the physical and functional multilayer network. 

``` r
ClinPriorScore<-MatrixPropagation(Y,alpha=0.2)

head(ClinPriorScore)
```


| Symbol |	geneID |	PriorFunct |	Symbol	| geneID | PriorPhys |
| --- | --- | --- | --- | --- | --- |
| L2HGDH |	79944 |	0.75624798 |	L2HGDH |	79944 |	0.762632978 |
| SUMF1 |	285362 |	0.673091841 |	SUMF1 |	285362 |	0.6754931 |
| NDUFAF4 |	29078 |	0.587444509 |	ADGRG1 |	9289 |	0.672186086 |
| NDUFA1 |	4694 |	0.578958636 |	NDUFAF4 |	29078 |	0.60198354 |
| ADGRG1 |	9289 |	0.548488778 |	NDUFA1 |	4694 |	0.586792265 |
| ASPA |	443 |	0.513182187 |	ASPA |	443 |	0.529834244 |

The best genes that match the HPO of the patient according to the propagation obtained in the functional (PriorFunct) and physical (PriorPhys) layers.
