# Analysis inputs

Place study inputs at the paths in `config/config.R` and `config/spatial_config.R`. Data are distributed separately from this code-only release. No large imaging files are included.

| Workflow | Required files |
|---|---|
| 01 | `supp_tables/Table_S4_segmented_PGP95_SYP_objects.xlsx`; `processed/clin_msi_matched.tsv`; `processed/clin_compact.tsv` |
| 02, 03, 06 | `processed/patient_nerve_quant.rds` (existing analysis-ready patient table) |
| 04 | `supp_tables/Table_S1_cohort_demographics.xlsx`; `external/BRS_cohort_A_counts.rds`; `external/BRS_cohort_B_counts.rds`; `external/BioMart_gene_anns.Rds`; the three GMT files named in config |
| 05 | `processed/scRNAseq_chen_annotated.rds`; `metadata/nerve_communication_genes_ann.rds`; DESeq2 results from 04 |
| 07 | Segmentation masks, phenotype assignments, reference measurement tables, image manifest and cohort annotations; see `docs/spatial/README.md` |
| 08 | `processed/spatial_measurements.rds`; `processed/patient_nerve_quant.rds`; the new Table S4 workbook |

RDS inputs use saveRDS/readRDS, except the BioMart archive, which contains the named workspace object `gene_anns`. Single-cell annotations require Gene, Class and Druggable columns; see `docs/SINGLE_CELL_HEATMAPS.md`.

The patient RDS remains necessary for clinical fields, molecular subtypes and stored STaN strata not fully represented in Table S1. No values are silently inferred from the publication table. Upstream clinical matching tables are not replaced by the downstream patient table.

Archive locations and access instructions must accompany the final manuscript data deposit. See the root README for existing data citations and current execution limits.
