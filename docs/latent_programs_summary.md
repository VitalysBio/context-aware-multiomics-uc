# Latent Biological Programs Identified in TCGA-BLCA

## Overview

Principal Component Analysis (PCA) was applied to the normalized RNA-seq expression matrix to derive latent transcriptomic programs capturing major sources of biological variation across bladder cancer patients.

Subsequent pathway enrichment, cluster association, mutation integration, CNV analysis, and survival modeling revealed that several latent components represent biologically interpretable tumor states with prognostic relevance.

These latent programs form the biological backbone of the Context-Aware Multi-Omics Patient Stratification framework.

---

## Latent Program Summary

| Program | Biological Interpretation                                               | Dominant Pathways                                                                                     | Key Marker Genes                                 | Associated Cluster(s)                     | Prognostic Role                    | Evidence                                                                                                       |
| ------- | ----------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------- | ------------------------------------------------ | ----------------------------------------- | ---------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| PC1     | Immune-stromal infiltration / macrophage-rich tumor microenvironment    | Immune System, Innate Immunity, Inflammatory Response, Complement, EMT, ECM Organization              | CD163, CSF1R, FCGR3A, FCGR2A, CD14, C1QC, FCER1G | Cluster 3 (highest), Cluster 2, Cluster 0 | Important survival feature in RSF  | Strong enrichment of myeloid and stromal pathways. Likely captures immune-excluded and macrophage-rich tumors. |
| PC3     | Proliferation and cell-cycle activity                                   | Cell Cycle, Mitotic Spindle, DNA Replication                                                          | Cell-cycle associated genes                      | Predominantly Cluster 2                   | Risk-associated                    | High proliferation signature commonly associated with aggressive tumors.                                       |
| PC4     | Cell motility and cytoskeletal remodeling                               | Migration, Cytoskeleton Organization, Cell Adhesion                                                   | Motility-related genes                           | Cluster 0 and Cluster 3                   | Risk-associated                    | Represents invasive cellular behavior and structural remodeling.                                               |
| PC6     | IFNγ-driven immune activation versus stromal EMT axis                   | Interferon Gamma Response, Interferon Alpha Response, Cytokine Signaling, Adaptive Immunity, EMT, ECM | Interferon-response genes                        | Cluster 3 and Cluster 0                   | High predictive importance         | Captures balance between anti-tumor immune activity and stromal remodeling.                                    |
| PC8     | RTK signaling, RHO GTPase activation, invasion and metastatic potential | Signal Transduction, RHO GTPase Cycle, RAC1 Signaling, Receptor Tyrosine Kinases, MET signaling       | Motility and signaling genes                     | Cluster 3                                 | Significant Cox predictor (HR > 1) | Represents invasive signaling programs linked to tumor progression.                                            |
| PC9     | Differentiation and epithelial organization                             | Cilium Assembly, Cilium Organization, Axonogenesis, Epithelial Differentiation                        | Differentiation-associated genes                 | Cluster 1 (highest), Cluster 3            | Protective (HR < 1)                | Consistent with differentiated luminal-like tumors and favorable prognosis.                                    |

---

## Relationship Between Latent Programs and Transcriptomic Clusters

### Cluster 0: Stromal / ECM-rich

Characteristics:

* EMT activation
* Extracellular matrix remodeling
* Collagen organization
* Intermediate PC1 levels
* Elevated stromal biology

---

### Cluster 1: Luminal-like Metabolic

Characteristics:

* Fatty acid metabolism
* Lipid metabolism
* Estrogen response
* Strongly negative PC1 values
* High epithelial differentiation

Most representative latent program:

* PC9

---

### Cluster 2: Basal / Squamous Inflammatory

Characteristics:

* TNF-alpha signaling
* Inflammatory response
* Keratinization
* Proliferative phenotype

Most representative latent program:

* PC3

---

### Cluster 3: Immune-Microenvironment Rich

Characteristics:

* Macrophage infiltration
* Complement activation
* Innate immune response
* EMT and ECM activity
* Highest PC1 values

Most representative latent programs:

* PC1
* PC6
* PC8

---

## Prognostic Programs Identified by Survival Modeling

The most reproducible prognostic signals across Cox and Random Survival Forest models were:

1. PC9 (Differentiation program)
2. PC6 (Immune activation program)
3. PC3 (Proliferation program)
4. PC1 (Immune-stromal microenvironment program)
5. Mutation burden
6. PC8 (Invasive signaling program)

These components consistently outperformed individual driver mutations and individual CNV events as survival predictors.

---

## Biological Interpretation

The major axes of transcriptomic variation in TCGA-BLCA appear to be driven by:

1. Tumor differentiation state.
2. Immune microenvironment composition.
3. Proliferative activity.
4. Stromal remodeling and EMT.
5. Invasive receptor tyrosine kinase signaling.

Together, these latent programs provide a biologically interpretable representation of bladder cancer heterogeneity and form the basis of the Context-Aware Multi-Omics Patient Stratification framework.
