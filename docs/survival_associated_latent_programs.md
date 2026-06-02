# Survival-associated latent transcriptomic programs

## Context

The integrated multi-omics Cox model identified PC8, PC9, and mutation burden as significant predictors of overall survival.

The full multi-omics Cox model achieved a C-index of 0.75, outperforming clinical-only, RNA-only, and mutation-only models.

## PC8: risk-associated motility and signaling program

PC8 was associated with increased hazard in the full multi-omics Cox model.

Main enriched pathways among PC8-positive genes included:

- RHO GTPase Cycle
- RAC1 GTPase Cycle
- RHOA GTPase Cycle
- Receptor Tyrosine Kinase Signaling
- MET Promotes Cell Motility
- Mitotic Spindle

Interpretation:

PC8 appears to capture a risk-associated program involving cytoskeletal remodeling, receptor tyrosine kinase signaling, RHO/RAC activity, and tumor cell motility.

PC8 differed significantly across RNA-seq clusters:

Kruskal-Wallis p = 1.19e-04

PC8 was highest in the immune-inflamed / EMT-associated cluster, suggesting that this latent risk axis may reflect an invasive and motility-associated state that cuts across discrete transcriptomic clusters.

## PC8-negative biology

PC8-negative genes were enriched for:

Mitochondrial respiratory chain complex assembly
Aerobic electron transport chain
Mitochondrial ATP synthesis coupled electron transport

### Interpretation:

Low PC8 may reflect a more oxidative metabolic state, whereas high PC8 may reflect a shift toward invasive signaling and motility programs.

## PC9: protective epithelial/ciliary organization program

PC9 was associated with reduced hazard in the full multi-omics Cox model.

Main enriched pathways among PC9-positive genes included:

Cilium Assembly
Cilium Organization
Anchoring of Basal Body to Plasma Membrane
Cell Projection Assembly
Axonogenesis / projection guidance

### Interpretation:

PC9 appears to capture a protective program related to epithelial organization, ciliary biology, and structured cellular projection organization.

PC9 differed significantly across RNA-seq clusters:

Kruskal-Wallis p = 1.63e-04

PC9 was lowest in the basal / squamous inflammatory cluster, suggesting loss of this organized epithelial/ciliary program may contribute to the aggressive behavior of this subtype.

## PC9-negative biology

PC9-negative genes were enriched for:

Myogenesis
Muscle contraction
KRAS signaling down
Estrogen response
Xenobiotic metabolism
mTORC1 signaling
Coagulation
Cholesterol homeostasis

### Interpretation:

The negative side of PC9 may represent stromal, smooth muscle, vascular, or metabolic microenvironmental signals rather than a purely tumor-intrinsic epithelial program.

Integrated interpretation

The survival-associated PCs suggest that prognosis in TCGA-BLCA is influenced by continuous latent transcriptomic programs in addition to discrete clusters.

PC8 represents a risk-associated motility and RTK/RHO signaling axis.

PC9 represents a protective epithelial/ciliary organization axis.

Together, these latent programs help explain why multimodal survival prediction improves when RNA-seq embeddings are integrated with mutation, CNV, and clinical features.

Working conceptual model
Cluster 0: stromal / ECM-rich, poor prognosis, likely driven by ECM and microenvironment remodeling.
Cluster 2: basal / squamous inflammatory, poor prognosis, low PC9 protective program, high CDKN2A deletion.
Cluster 3: immune-inflamed / EMT-high, high PC8 but also high PC9, suggesting mixed invasive and organized/ciliary biology.
Cluster 1: luminal-like metabolic, relatively better prognosis and intermediate/protective latent features.

## Next step

The next modeling step will compare penalized Cox regression with non-linear survival models such as Random Survival Forest, followed by feature importance and uncertainty-aware prediction.