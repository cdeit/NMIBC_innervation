# Reproduction status

## Validation performed

Clinical (02), subtype plots (03), bulk RNA-seq (04) and nerve-cluster analyses (06) completed with the analysis-ready inputs in the development environment. The bulk run used 97 RNA-matched patients, selecting 33 STaN-low and 33 STaN-high samples; the selected IDs match the combined gotNeRve training/test sample IDs. Script 04 now reads Table S1: its IDs, covariates and abundance values were compared with the earlier RDS input. Models were not refitted for that input-routing change.

The adequate-BCG subset of the expanded patient RDS reproduces the earlier 141-patient table. Spatial models and bubble plots previously completed on archived measurements. S4F mapping preserves the same 45 eligible spatial images and stored patient strata. These checks do not establish equivalence of every figure or statistic to the manuscript.

## Remaining requirements

| Workflow | Remaining requirement |
|---|---|
| 01 | Historical matched/compact clinical tables; the complete S4 workbook is supported. Full phenotyping/aggregation has not been rerun from these publication inputs. |
| 02, 03, 06 | Analysis-ready patient RDS remains required for fields absent from S1. Optional imaging-survival mode in 03 is not configured. |
| 04 | S1, A/B counts, BioMart archive and three GMT files; exact input paths are in config/config.R. |
| 05 heatmaps | Final annotated Chen object and nerve_communication_genes_ann.rds. Syntax and annotation-reader checks passed; full heatmaps have not been regenerated. Missing-gene/column-split alignment and cell-type completeness require validation with the reference object. |
| 05 preprocessing | Reconcile the cleaned multiplet-rate table with provenance/process_singlecell_original.R before execution. The cleaned table currently pairs 11 rates with 20 cell-count entries. The optional preprocessing command stops explicitly. |
| 07 | External TIFFs, phenotype assignments, image manifest, clinical annotations, reference tables and preserved Alpine environment. A five-image end-to-end pass was reported; full-cohort reproduction has not been established here. |
| 08 | Archived measurements, patient RDS and S4; distance-distribution figures additionally require the original distribution_palette in config/spatial_config.R. Models can complete before this plotting stage stops. |

Table S1 does not fully replace the patient RDS. In particular, its nine LVI entries labelled NA; TaHG correspond to three No and six missing values in the earlier table. Stored STaN groups and some subtype fields are absent. Existing values remain in the RDS rather than being reconstructed from ambiguous labels.

## Environment

The recorded downstream runs used R 4.3.3; recorded versions include DESeq2 1.42.1, clusterProfiler 4.10.1, GOSemSim 2.28.1, rstatix 1.1.0, dplyr 1.2.1 and ggplot2 4.0.3. This is a partial environment record, not a complete lockfile. Library calls and namespace-qualified dependencies are explicit in the code. A complete reference environment is still needed for independent reproduction.

Spatial image processing uses a separate R 4.4.1/Linux container and numerical-library requirements; see [spatial system requirements](spatial/SYSTEM_REQUIREMENTS.md). Do not use downstream plotting package versions as a substitute for that runtime. No dependencies are installed or upgraded by the scripts.

## Repository checks

The release is checked for R syntax, shell syntax, configuration/source paths and local documentation links. Those checks do not execute scientific analyses. Input archives, the exact environment and reference outputs are required for full reproduction.
