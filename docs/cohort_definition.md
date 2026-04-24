# Cohort Definition

## Purpose

This document defines the analytical cohort for the project and establishes the unit of analysis, inclusion logic, and initial modality scope.

## 1. Disease Context

- Disease: urothelial carcinoma / bladder cancer
- Discovery cohort: TCGA-BLCA
- External validation cohort: IMvigor210

## 2. Unit of Analysis

The primary unit of analysis will be the **patient**.

When multiple samples exist for the same patient, selection rules will be defined to retain a single representative tumor sample for modeling unless a later analysis explicitly addresses multi-sample structure.


## 3. Initial Discovery Cohort Scope

The first implementation will focus on TCGA-BLCA patients with at least:

- clinical metadata
- RNA-seq data

Additional modalities to be integrated when available:

- somatic mutations
- copy number variation

This strategy allows the project to start with a strong transcriptomic backbone while preserving multi-omics expansion.


## 4. Initial Inclusion Criteria

Patients will be eligible for the initial analytical cohort if they meet all of the following:

- valid TCGA patient identifier
- available clinical metadata
- available RNA-seq data
- tumor sample suitable for analysis

For multimodal analyses, patients will additionally require the corresponding modality-specific data.

## 5. Initial Exclusion Criteria

Patients or samples may be excluded if they have:

- invalid or non-matching identifiers
- duplicated records
- non-primary tumor samples where primary tumor is preferred
- missing essential metadata
- failed quality control

## 6. Sample Selection Principles

When multiple biospecimen entries exist, priority will generally be given to:

1. primary tumor samples
2. samples with complete modality coverage
3. samples with the clearest metadata consistency

Exact rules will be finalized during data harmonization.

## 7. Identifier Harmonization

Expected identifier levels include:

- file-level identifiers
- sample-level identifiers
- aliquot-level identifiers
- TCGA barcode-derived patient identifiers

A harmonization layer will be required to map all files to a patient-level master table.


## 8. Initial Modeling Cohorts

### Cohort A: Discovery / Stratification
TCGA-BLCA patients with clinical + RNA-seq  
Use:
- exploratory analysis
- embeddings
- clustering
- survival-linked stratification

### Cohort B: Expanded Multi-Omics Discovery
TCGA-BLCA patients with clinical + RNA-seq + mutation + CNV  
Use:
- multimodal modeling
- ablation analysis
- context-aware fusion development

### Cohort C: External Therapeutic Validation
IMvigor210 patients with expression + clinical + response labels  
Use:
- external validation for immunotherapy response prediction
- uncertainty and calibration assessment in a clinically relevant setting


## 9. Rationale for This Cohort Design

This staged cohort strategy balances:

- scientific rigor
- practical feasibility
- robustness to missing modalities
- portfolio value
- external clinical relevance

It avoids delaying the project until every modality is perfectly harmonized while preserving a clear roadmap toward true multimodal analysis.
