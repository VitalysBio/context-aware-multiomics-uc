# Scientific Rationale

## 1. Background

Urothelial carcinoma (bladder cancer) is characterized by high molecular heterogeneity across patients, involving alterations at multiple biological levels, including gene expression, somatic mutations, copy number variation, and clinical features.

Immune checkpoint inhibitors have improved outcomes in a subset of patients; however, response rates remain limited and highly variable. Current biomarkers such as PD-L1 expression or tumor mutational burden (TMB) provide only partial predictive power.

This suggests that treatment response is not driven by a single modality but rather by complex, patient-specific interactions across multiple biological layers.

## 2. Problem Statement

Most existing computational approaches for multi-omics integration rely on:
- early concatenation of features
- modality-agnostic embeddings
- static weighting of data sources

These approaches assume that all modalities contribute equally to prediction across all patients, which is biologically unrealistic.

Additionally, many models provide point predictions without any estimate of uncertainty, limiting their usefulness in clinical decision-making.


## 3. Central Hypothesis

Treatment response in urothelial carcinoma emerges from heterogeneous combinations of molecular and clinical signals, and these combinations vary across patients.

A model that can dynamically adjust the contribution of each modality at the patient level will better capture this heterogeneity.


## 4. Modeling Hypothesis

A context-aware multi-omics model that learns patient-specific modality weights will:
- outperform static integration approaches
- improve calibration of predictions
- provide interpretable insights into modality relevance


## 5. Role of Uncertainty

In a clinical context, it is not sufficient to produce accurate predictions. It is equally important to quantify confidence.

We hypothesize that:
- uncertain predictions will correspond to biologically ambiguous or out-of-distribution cases
- uncertainty-aware modeling will improve reliability and interpretability


## 6. Biological Interpretation Goals

This project aims to go beyond feature importance and enable:

- identification of modality-specific drivers of prediction
- discovery of biologically meaningful patient subgroups
- pathway-level interpretation of predictive signals
- comparison of molecular programs across responders and non-responders


## 7. Translational Relevance

The project is designed to align with real-world applications in precision medicine, including:

- patient stratification
- treatment selection
- biomarker discovery
- clinical decision support systems

## 8. Differentiation from Standard Projects

This project intentionally avoids:

- simple concatenation-based integration
- purely descriptive clustering without clinical linkage
- black-box modeling without interpretability
- evaluation without uncertainty or calibration

Instead, it focuses on a structured, interpretable, and clinically grounded multi-omics modeling framework.
