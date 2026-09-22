# Expression context for bulk STaN-associated and curated neural-TME genes.
# Chen samples are not assigned STaN-high/low labels.

#' Map bulk STaN-associated genes onto annotated NMIBC cell types
#' @param sc_file Annotated Chen Seurat object, stored with saveRDS().
#' @param deg_file Bulk DESeq2 result CSV.
#' @param annotation_file Prepared saveRDS() table with Gene, Class and Druggable columns.
#'   Row order defines the curated gene order; annotations are not derived here.
#' @param out_dir Figure output directory.
run_single_cell_heatmaps <- function(sc_file, deg_file, annotation_file,
                                    out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  # ---- Load NMIBC single-cell reference ----------------------------------------
  sc <- readRDS(sc_file)
  DefaultAssay(sc) <- "RNA"
  Idents(sc) <- "cell_type"

  sc_nmibc <- subset(sc, subset = Tumor_type == "NMIBC" & Tumor_stage == "NMIBC")
  message(
    "NMIBC scRNA-seq reference: ",
    length(unique(sc_nmibc$Seq_ID)), " SRA accessions; ", ncol(sc_nmibc), " cells."
  )

  cell_order <- c(
    "Urothelial cell", "Proliferating cell", "iCAF", "MyoCAF", "Endothelial",
    "Myeloid", "Differentiating plasma cell", "Mature plasma cell",
    "B cell", "T cell", "Mast cell"
  )

  # ---- Shared heatmap function -------------------------------------------------
  # Average normalized RNA expression is calculated for each annotated cell type.
  # Genes are columns and cell types are rows.
  plot_celltype_heatmap <- function(sc, genes, annotations = NULL, column_split = NULL,
                                    title = NULL, cluster_within_slices = TRUE) {
    genes <- unique(genes[genes %in% rownames(sc)])
    if (!length(genes)) stop("None of the requested genes are present in the Seurat object.")

    avg <- AverageExpression(
      sc,
      features = genes, assays = "RNA", layer = "data",
      group.by = "cell_type", verbose = FALSE
    )$RNA

    present_cell_types <- cell_order[cell_order %in% colnames(avg)]
    avg <- t(as.matrix(avg[genes, present_cell_types, drop = FALSE]))

    col_fun <- colorRamp2(
      c(0, max(avg, na.rm = TRUE)),
      hcl_palette = "Purples 3",
      reverse = TRUE
    )

    top_ann <- NULL
    if (!is.null(annotations)) {
      annotations <- annotations[match(colnames(avg), annotations$Gene), , drop = FALSE]
      if (anyNA(annotations$Gene))
        stop("Annotations are missing for one or more plotted genes.")

      ann_data <- annotations[, setdiff(names(annotations), "Gene"), drop = FALSE]
      top_ann <- HeatmapAnnotation(df = ann_data)
    }

    Heatmap(
      avg,
      name = "Average\nExpression",
      col = col_fun,
      cluster_rows = FALSE,
      row_order = present_cell_types,
      show_row_names = TRUE,
      show_row_dend = FALSE,
      cluster_columns = FALSE,
      column_split = column_split,
      cluster_column_slices = cluster_within_slices,
      clustering_method_columns = "ward.D2",
      show_column_names = TRUE,
      show_column_dend = FALSE,
      column_names_side = "bottom",
      column_names_rot = 45,
      top_annotation = top_ann,
      column_title = title,
      border = TRUE
    )
  }

  save_heatmap <- function(ht, filename, width, height = 7) {
    pdf(filename, width = width, height = height)
    draw(
      ht,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      merge_legends = TRUE
    )
    dev.off()
    message("Saved: ", filename)
  }

  # ==============================================================================
  # 1. BULK STaN-ASSOCIATED DEG HEATMAP
  # ==============================================================================

  degs <- read.csv(deg_file, check.names = FALSE)

  # These are the significant bulk STaN-high vs STaN-low DEGs used in the paper
  # (FDR < 0.05 and |log2FC| > 1). The single-cell samples are not STaN-stratified.
  deg_plot <- degs |>
    filter(
      !is.na(padj),
      padj < 0.05,
      abs(log2FoldChange) > 1
    ) |>
    mutate(
      Direction = if_else(
        log2FoldChange > 0,
        "Upregulated in STaN-high",
        "Downregulated in STaN-high"
      )
    ) |>
    arrange(desc(log2FoldChange))

  deg_genes <- deg_plot$Gene
  deg_split <- factor(
    deg_plot$Direction,
    levels = c("Upregulated in STaN-high", "Downregulated in STaN-high")
  )

  ht_deg <- plot_celltype_heatmap(
    sc_nmibc,
    genes = deg_genes,
    column_split = deg_split,
    title = "Bulk STaN-associated DEGs"
  )

  save_heatmap(
    ht_deg,
    file.path(out_dir, "scRNAseq_celltype_heatmap_STaN_DEGs.pdf"),
    width = max(12, length(deg_genes) * 0.12)
  )

  # ==============================================================================
  # 2. CURATED NEURAL-TME COMMUNICATION HEATMAP
  # ==============================================================================

  # Read final custom/Human Protein Atlas annotations without reclassifying genes.
  # One row per curated gene; preserve the supplied order and annotation labels.
  curated_ann <- readRDS(annotation_file)
  required <- c("Gene", "Class", "Druggable")
  if (!is.data.frame(curated_ann) || !all(required %in% names(curated_ann))) {
    stop("Annotation RDS must contain a data frame with columns: Gene, Class, Druggable.")
  }
  curated_ann <- curated_ann[, required, drop = FALSE]
  if (!nrow(curated_ann) || anyNA(curated_ann) ||
      any(vapply(curated_ann, function(x) any(!nzchar(x)), logical(1))) ||
      anyDuplicated(curated_ann$Gene)) {
    stop("Provide one complete annotation row per gene; duplicate/blank values are not allowed.")
  }
  curated_genes <- curated_ann$Gene

  # Preserve the original annotation palettes.
  class_cols <- c(
    "Adrenergic receptor" = "deeppink",
    "Neuropeptide receptor" = "#C4B8F2",
    "Glutamatergic receptor" = "mediumpurple",
    "Neurotrophic factor" = "orange",
    "Neural adhesion molecule" = "tan",
    "Postsynaptic scaffold/signaling" = "deepskyblue",
    "Axon growth/attraction" = "lightgreen",
    "Axon repulsion" = "mediumvioletred"
  )
  drug_cols <- c(
    "Potential drug target" = "gold",
    "FDA approved drug target" = "mediumaquamarine",
    "No/Unknown" = "lightgrey"
  )

  if (!all(curated_ann$Class %in% names(class_cols)) ||
      !all(curated_ann$Druggable %in% names(drug_cols))) {
    stop("Annotation labels must match the established palettes; see docs/SINGLE_CELL_HEATMAPS.md.")
  }

  # Build this heatmap explicitly so the final figure retains its annotation colors.
  genes_present <- curated_genes[curated_genes %in% rownames(sc_nmibc)]
  avg_curated <- AverageExpression(
    sc_nmibc,
    features = genes_present, assays = "RNA", layer = "data",
    group.by = "cell_type", verbose = FALSE
  )$RNA
  avg_curated <- t(as.matrix(avg_curated[genes_present, cell_order, drop = FALSE]))

  ann_curated <- curated_ann[match(colnames(avg_curated), curated_ann$Gene), ]
  curated_split <- factor(ann_curated$Class)

  top_ann <- HeatmapAnnotation(
    Druggable = ann_curated$Druggable,
    Class = ann_curated$Class,
    height = unit(2, "mm"),
    col = list(Druggable = drug_cols, Class = class_cols),
    annotation_legend_param = list(
      Druggable = list(title = "Druggable", position = "bottom"),
      Class = list(title = "Class", position = "bottom")
    )
  )

  curated_col_fun <- colorRamp2(
    c(0, max(avg_curated, na.rm = TRUE)),
    hcl_palette = "Purples 3",
    reverse = TRUE
  )

  ht_curated <- Heatmap(
    avg_curated,
    name = "Average\nExpression",
    col = curated_col_fun,
    cluster_rows = FALSE,
    row_order = cell_order,
    show_row_names = TRUE,
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    column_split = curated_split,
    cluster_column_slices = TRUE,
    clustering_method_columns = "ward.D2",
    show_column_names = TRUE,
    show_column_dend = FALSE,
    column_names_side = "bottom",
    column_names_rot = 45,
    top_annotation = top_ann,
    column_title = "Candidate neural-TME communication genes",
    border = TRUE
  )

  save_heatmap(
    ht_curated,
    file.path(out_dir, "fig7_scRNAseq_neural_TME_heatmap.pdf"),
    width = max(12, length(genes_present) * 0.12)
  )
  invisible(list(deg_heatmap = ht_deg, curated_heatmap = ht_curated))
}
