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

5.77e-03

This suggests that transcriptomic stratification captures clinically relevant tumor biology.

Cluster 0: Stromal / ECM-rich / EMT-like

Representative marker genes
MOXD1
PODN
CDH11
EMILIN1
COL14A1
LUM
SMOC2
DCN
COL8A1
FGFR1

Enriched pathways

Epithelial Mesenchymal Transition
Extracellular Matrix Organization
Collagen Biosynthesis
Collagen Chain Trimerization
Collagen Fibril Organization
Regulation of Angiogenesis
Interpretation

Cluster 0 appears to represent a stromal-rich and extracellular matrix-associated tumor state. The strong enrichment of collagen organization, EMT, and angiogenesis pathways suggests a microenvironmental or fibroblast-associated program.

Working label:

Stromal / ECM-rich / EMT-like

Cluster 1: Luminal-like metabolic / differentiated epithelial

Representative marker genes
TBX3
RNF128
GATA3-AS1
VSIG2
DHRS2
SRCIN1
RAB15

Enriched pathways

Metabolism of Lipids
Fatty Acid Metabolism
Arachidonic Acid Metabolism
Estrogen Response Early / Late
Epithelium Development
Epithelial Cell Differentiation
Xenobiotic Metabolism
Cytochrome P450-related pathways

Interpretation

Cluster 1 appears compatible with a luminal-like or differentiated epithelial phenotype. The enrichment of lipid metabolism, epithelial differentiation, estrogen response, and xenobiotic metabolism suggests a more metabolically differentiated tumor state.

Working label:

Luminal-like metabolic

Cluster 2: Basal / squamous inflammatory

Representative marker genes
KRT6A
KRT5
DSP
SERPINB5
DSG3
CD44
IL1RAP
TFAP2A

Enriched pathways

Epidermis Development
Keratinization
Formation of Cornified Envelope
TNF-alpha Signaling via NF-kB
Inflammatory Response
Interferon Gamma Response
Hypoxia
Epithelial Mesenchymal Transition

Interpretation

Cluster 2 shows a strong basal/squamous-like phenotype supported by keratinization, epidermal development, and basal epithelial marker genes. The additional enrichment of inflammatory, hypoxia, and EMT pathways suggests a potentially aggressive and plastic tumor state.

Working label:

Basal / squamous inflammatory

Cluster 3: Immune-inflamed / IFNγ-high / EMT-associated

Representative marker genes
THEMIS2
SIRPA
ITGB2
PIK3R5
WIPF1
SLC7A7
ANXA6

Enriched pathways
Allograft Rejection
Immune System
Interferon Gamma Response
Inflammatory Response
Cytokine Signaling in Immune System
Innate Immune System
Adaptive Immune System
Regulation of T Cell Activation
Epithelial Mesenchymal Transition
Extracellular Matrix Organization

Interpretation

Cluster 3 appears to represent an immune-inflamed tumor state with strong IFNγ, cytokine signaling, innate/adaptive immune activation, and T-cell-related signatures. The concurrent EMT and extracellular matrix enrichment suggests an immune-inflamed but mesenchymal-associated microenvironment.

Working label:

Immune-inflamed IFNγ-high

Interpretation caveats

These labels are preliminary and based on transcriptomic clustering, marker gene analysis, pathway enrichment, and survival association. They should be interpreted as biologically informed working labels rather than definitive molecular subtype calls.

Future work should validate these assignments using:

published BLCA subtype signatures
mutation and CNV integration
immune/stromal scoring
external immunotherapy cohorts
