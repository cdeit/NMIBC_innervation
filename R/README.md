# Analysis functions

Run numbered scripts from the repository root; do not source all function files to launch analyses. The [root README](../README.md) maps each workflow to its inputs and supporting files.

- `read_table_s4.R`: supplementary workbook reader with required-sheet checks.
- `nerve_clustering.R`: Stage 1 identification, Stage 2 published/de novo phenotyping and diagnostic/manuscript embeddings.
- `nerve_aggregation.R`: object, region, image and patient filtering and aggregation.
- `clinical_plots.R`, `nerve_cluster_analysis.R`, `survival.R`: clinical/subtype/cluster plots and survival analyses.
- `DESeq2_and_GSEA.R`: bulk DESeq2 and gene-set enrichment.
- `singlecell_preprocessing.R`, `singlecell_heatmaps.R`: Chen processing and bulk-DEG cell-type context, preprocessing and heatmaps in separate functions.
- `prepare_spatial_masks.R`, `compute_spatial_measurements.R`: active segmentation preparation and spatial distance/band measurements, orchestrated by `spatial_image_workflow.R`.
- `spatial_analysis.R`: downstream spatial models and plots.
- `data_io.R`: shared input-file checks and loading named objects from supplied workspace archives.

General settings are in `config/config.R`; spatial settings remain in `config/spatial_config.R`. Sourcing either config does not create output directories. Existing external identifiers and quantitative column names are retained.

Published Stage 2 labels are selected by `cfg$clustering$use_publication_cluster_labels`. De novo runs are reanalyses, not a promise of exact manuscript label reproduction. Scaling remains based on the broader Stage 1 population. Image eligibility is applied during aggregation, not moved upstream.

See [review items](../docs/REPRODUCIBILITY.md) for missing inputs, ambiguous historical code and known execution blockers. Do not infer scientific replacements from these issues.
