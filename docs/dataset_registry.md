# Dataset Registry

## Purpose

This document tracks all datasets used in the project, including their source, modality coverage, intended use, and known limitations.

This ensures:
- reproducibility
- traceability of decisions
- clarity of data provenance

## 1. TCGA-BLCA (Discovery Cohort)

### Source
- The Cancer Genome Atlas (TCGA)
- Accessed via GDC portal

### Modalities
- RNA-seq (gene expression)
- somatic mutations (MAF files)
- copy number variation (CNV)
- clinical metadata
- survival data

### Intended Use
- multi-omics integration
- embedding learning
- clustering and stratification
- survival analysis
- feature engineering

### Strengths
- large cohort size
- multi-omics coverage
- standardized preprocessing
- rich clinical annotations

### Limitations
- no direct immunotherapy response labels
- batch effects across centers
- missing data in some modalities


## 2. IMvigor210 (External Validation Cohort)

### Source
- Public immunotherapy study (atezolizumab-treated patients)
- Available via cBioPortal or associated R packages

### Modalities
- gene expression (primarily RNA-based)
- clinical data
- treatment response labels

### Intended Use
- evaluation of predictive performance for immunotherapy response
- external validation of model generalization
- clinical relevance assessment

### Strengths
- real treatment response data
- clinically meaningful endpoint
- widely used benchmark in immunotherapy studies

### Limitations
- limited multi-omics compared to TCGA
- potential differences in preprocessing pipelines
- smaller sample size


## 3. Gene Set Databases (Later Phase)

### Sources
- MSigDB
- Reactome
- Hallmark gene sets

### Intended Use
- pathway aggregation
- enrichment analysis
- biological interpretation


## 4. Inclusion Criteria (Initial Plan)

Patients will be included if:
- they have valid identifiers across modalities
- they have at least one key modality (RNA-seq preferred)
- clinical metadata is available

Further filtering rules will be defined during preprocessing.


## 5. Exclusion Criteria (Initial Plan)

- samples with missing or inconsistent IDs
- duplicated entries
- samples without usable expression data
- samples failing quality control


## 6. Data Harmonization Challenges

Expected issues:
- inconsistent patient identifiers across modalities
- gene ID formats (Ensembl vs gene symbols)
- missing modalities per patient
- batch effects

Strategies will be defined in the preprocessing phase.

## 7. Versioning Strategy

Each dataset will be:
- stored in raw, interim, and processed formats
- associated with configuration files
- documented with download scripts


## 8. Future Dataset Extensions (Optional)

Potential additions:
- CPTAC proteomics
- GEO external cohorts
- additional immunotherapy datasets

These will only be added if they provide clear value and do not introduce excessive complexity early in the project.
