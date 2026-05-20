# Somatic Mutation Integration Notes

## Objective

Integrate somatic mutation data from TCGA-BLCA to characterize the genomic landscape underlying RNA-seq-derived transcriptomic tumor states.

## Data source

Somatic mutation data were downloaded from the GDC portal using:

- Project: TCGA-BLCA
- Data Category: Simple Nucleotide Variation
- Data Type: Masked Somatic Mutation
- Workflow Type: Aliquot Ensemble Somatic Variant Merging and Masking
- Access: Open

## Processing strategy

Mutation Annotation Format files were parsed and filtered to retain non-silent coding variants, including missense, nonsense, frameshift, splice-site, in-frame indels, nonstop mutations, and translation start site mutations.

Silent and non-coding/modifier variants were excluded to focus on functionally relevant alterations.

## Outputs generated

- `tcga_blca_mutation_table_long.csv`
- `tcga_blca_mutation_burden.csv`
- `tcga_blca_mutation_binary_matrix.csv`
- `recurrent_driver_genes.csv`
- `driver_mutation_enrichment_by_cluster.csv`

## Main results

A total of 84,746 non-silent mutations were detected across 406 patients.

The most recurrently mutated genes included:

- TP53
- KMT2D
- KDM6A
- ARID1A
- PIK3CA
- RB1
- EP300
- KMT2C

Mutation burden did not differ significantly across RNA-seq clusters:


Kruskal-Wallis p = 0.218

This suggests that transcriptomic tumor states are not simply driven by global mutation burden.

## Driver mutation patterns by cluster

Cluster 1, labeled as luminal-like metabolic, showed relative depletion of TP53 mutations.

Cluster 2, labeled as basal / squamous inflammatory, showed high TP53 mutation frequency.

Cluster 3, labeled as immune-inflamed IFNγ-high, showed high TP53 frequency, RB1 enrichment, and KDM6A depletion.

## Interpretation

The mutation layer supports the biological relevance of the RNA-seq-derived clusters. The results suggest that transcriptomic tumor states are associated with distinct driver mutation landscapes, rather than reflecting only differences in overall mutation burden.

## Next step

The next planned omics layer is copy number variation integration, followed by context-aware multi-omics fusion.

