# Preliminary RNA-seq Cluster Interpretation Notes

## Context

RNA-seq based clustering of TCGA-BLCA primary tumor samples identified four transcriptomic clusters.

The clusters were generated using:
- log2 CPM-normalized RNA-seq expression
- highly variable genes
- PCA dimensionality reduction
- KMeans clustering

Kaplan-Meier analysis showed significant survival differences across clusters using a global log-rank test.

## Survival association

Global log-rank p-value:

```text
5.77e-03```

This suggests that transcriptomic stratification captures clinically relevant tumor biology.

Cluster 0: Stromal / extracellular matrix-rich phenotype

Top markers include:

MOXD1
PODN
CDH11
EMILIN1
COL14A1
LUM
SMOC2
TMEM119

Preliminary interpretation:

This cluster appears enriched for extracellular matrix, collagen organization, fibroblast-like, and stromal remodeling signals.

Working label:
Stromal / ECM-rich

Cluster 1: Luminal-like / differentiated urothelial phenotype

Top markers include:

TBX3
RNF128
GATA3-AS1
VSIG2
DHRS2
SRCIN1

Preliminary interpretation:

This cluster appears compatible with a more luminal-like or differentiated urothelial phenotype.

Working label:

Luminal-like

Cluster 2: Basal / squamous-like phenotype

Top markers include:

KRT6A
KRT5
DSP
SERPINB5
DSG3
CD44
TFAP2A
IL1RAP

Preliminary interpretation:

This cluster shows a strong basal/squamous-like program, supported by keratin markers and epithelial basal markers.

Working label:

Basal / squamous-like

Cluster 3: Immune / myeloid-enriched phenotype

Top markers include:

THEMIS2
SIRPA
ITGB2
PIK3R5
WIPF1
SLC7A7

Preliminary interpretation:

This cluster appears enriched for immune and potentially myeloid-associated signals.

Working label:

Immune / myeloid-enriched
