# Single-cell heatmaps

`R/singlecell_heatmaps.R` defines `run_single_cell_heatmaps()`, called by script 05. Preprocessing remains separate in `R/singlecell_preprocessing.R`; its execution settings still require validation.

Run from the repository root:

```sh
Rscript scripts/05_run_single_cell_analysis.R
```

Inputs:

- `cfg$paths$single_cell`: annotated Chen Seurat object, with RNA assay and `cell_type`, `Tumor_type`, `Tumor_stage`, `Seq_ID` metadata.
- Script 04's DESeq2 CSV under `results/bulk_RNAseq/tables/`.
- `cfg$paths$single_cell_annotations`: `data/metadata/nerve_communication_genes_ann.rds`, a prepared data frame saved with `saveRDS()` and read with `readRDS()`.

The annotation file has one row per curated gene and three required columns: `Gene`, `Class`, `Druggable`. Additional columns may be retained in the file but are not plotted. Row order supplies the curated gene order before the existing heatmap split/clustering settings. Gene symbols must be unique. Missing annotations are reported, not converted into negative classifications.

`data/metadata/single_cell_gene_annotations_template.tsv` preserves the previous 156 curated genes, their order and their Class assignments. Druggable is intentionally blank: fill it from the final annotations and save the completed data frame with `saveRDS(annotations, "data/metadata/nerve_communication_genes_ann.rds")`. The template is not a runnable final annotation file.

Established Class labels:

- Adrenergic receptor
- Neuropeptide receptor
- Glutamatergic receptor
- Neurotrophic factor
- Postsynaptic scaffold/signaling
- Axon growth/attraction
- Axon repulsion
- Neural adhesion molecule

Established Druggable labels: `Potential drug target`, `FDA approved drug target`, `No/Unknown`.

The plotter does not parse raw Human Protein Atlas fields, infer druggability, or rebuild curated classes. To preserve the existing figure, supply the same genes, order and annotation labels. A different annotation schema or additional tracks needs an explicit plotting update.

Outputs remain `results/single_cell/figures/scRNAseq_celltype_heatmap_STaN_DEGs.pdf` and `results/single_cell/figures/fig7_scRNAseq_neural_TME_heatmap.pdf`. The first heatmap uses bulk FDR < 0.05 and absolute log2 fold change > 1; the prepared Class/Druggable tracks apply to the second heatmap. Both use the existing NMIBC subset and average RNA-expression plotting code. Chen cells are not assigned STaN-high/low status.

Validation: parsed all modified R files; confirmed unchanged parsed expressions for DEG selection, heatmap construction, palettes, expression calculation and save helper; tested annotation ordering and duplicate rejection with a temporary synthetic fixture. Full figures have not been regenerated: the final annotation file and configured annotated Seurat input are not yet available. Existing missing-gene/split and cell-type completeness issues still require validation on the reference object.
