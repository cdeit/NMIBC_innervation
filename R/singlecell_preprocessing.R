# Chen et al. (2020), PRJNA662018: optional preprocessing of 10x matrices.
# STaN status is not measured in these samples. See docs/REPRODUCIBILITY.md for
# unresolved preprocessing settings; use the annotated reference for mapping.

#' Preprocess and annotate the Chen single-cell reference
#' @param data_root Directory containing SRR12603780-SRR12603790 10x outputs.
#' @param out_file Destination for the annotated Seurat RDS object.
run_single_cell_preprocessing <- function(data_root, out_file) {
  accessions <- paste0("SRR126037", 80:90)
  data_dirs <- file.path(data_root, accessions, "outs", "filtered_feature_bc_matrix")

  # Sample annotations obtained from the dataset authors, keyed directly by
  # SRA accession to avoid maintaining a redundant sample alias.
  sample_metadata <- data.frame(
    Seq_ID = c(
      "SRR12603790", "SRR12603789", "SRR12603787", "SRR12603786",
      "SRR12603785", "SRR12603784", "SRR12603783", "SRR12603782",
      "SRR12603781", "SRR12603780", "SRR12603788"
    ),
    Tumor_type = c(
      "Normal", "Normal", "MIBC", "NMIBC", "NMIBC", "MIBC",
      "NMIBC", "MIBC", "Normal", "MIBC", "NMIBC"
    ),
    Tumor_grade = c(
      "Adj. high grade", "Adj. high grade", "High grade", "High grade",
      "Low grade", "High grade", "High grade", "High grade",
      "Adj. low grade", "High grade", "Low grade"
    ),
    Tumor_stage = c(
      "MIBC", "MIBC", "MIBC", "NMIBC", "NMIBC", "MIBC",
      "NMIBC", "MIBC", "NMIBC", "MIBC", "NMIBC"
    ),
    stringsAsFactors = FALSE
  )

  # ---- Read and filter individual samples --------------------------------------
  message("Reading 10x matrices...")
  future::plan(multisession, workers = max(1, parallel::detectCores() - 1))

  seurat_list <- future.apply::future_lapply(seq_along(data_dirs), function(i) {
    x <- CreateSeuratObject(counts = Read10X(data.dir = data_dirs[i]))
    x$percent.mt <- PercentageFeatureSet(x, pattern = "^MT-")
    x <- subset(
      x,
      subset = nFeature_RNA > 300 & nFeature_RNA < 6000 &
        percent.mt < 10 & nCount_RNA > 1000
    )
    x$Seq_ID <- accessions[i]
    x
  })
  future::plan(sequential)

  post_filter_cellcount <- sum(vapply(seurat_list, ncol, integer(1)))
  message("Cells after expression/mitochondrial QC: ", post_filter_cellcount)

  # ---- Cell-cycle scoring and doublet detection --------------------------------
  s_genes <- cc.genes.updated.2019$s.genes
  g2m_genes <- cc.genes.updated.2019$g2m.genes |>
    gsub(pattern = "PIMREG", replacement = "FAM64A") |>
    gsub(pattern = "JPT1", replacement = "HN1")

  multiplet_rates_10x <- data.frame(
    Multiplet_rate = c(.004, .008, .016, .023, .031, .039, .046, .054, .061, .069, .076),
    Recovered_cells = seq(500, 10000, by = 500)
  )

  future::plan(multisession, workers = max(1, parallel::detectCores() - 1))
  options(future.globals.maxSize = 100 * 1024^3)

  seurat_list <- future.apply::future_lapply(seurat_list, function(x) {
    x <- NormalizeData(x)
    x <- ScaleData(x)
    x <- CellCycleScoring(x, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
    x <- FindVariableFeatures(x)
    x <- RunPCA(x, nfeatures.print = 10)

    # Select PCs using the criteria applied in the original analysis.
    pct <- x[["pca"]]@stdev / sum(x[["pca"]]@stdev) * 100
    co1 <- which(cumsum(pct) > 90 & pct < 5)[1]
    co2 <- sort(which((pct[-length(pct)] - pct[-1]) > 0.1), decreasing = TRUE)[1] + 1
    n_pcs <- min(co1, co2)

    x <- RunUMAP(x, dims = 1:n_pcs)
    x <- FindNeighbors(x, dims = 1:n_pcs)
    x <- FindClusters(x, resolution = 0.1)

    sweep_stats <- summarizeSweep(paramSweep(x, PCs = 1:n_pcs, sct = FALSE))
    bcmvn <- find.pK(sweep_stats)
    pK <- bcmvn |>
      filter(BCmetric == max(BCmetric)) |>
      slice(1) |>
      pull(pK) |>
      as.character() |>
      as.numeric()

    # Estimate the expected doublet rate from the 10x recovered-cell table.
    multiplet_rate <- multiplet_rates_10x |>
      filter(Recovered_cells < ncol(x)) |>
      slice_max(Recovered_cells, n = 1, with_ties = FALSE) |>
      pull(Multiplet_rate)

    homotypic_prop <- modelHomotypic(x$seurat_clusters)
    n_exp <- round(multiplet_rate * ncol(x) * (1 - homotypic_prop))

    x <- doubletFinder(x, PCs = 1:n_pcs, pK = pK, nExp = n_exp)
    df_col <- grep("^DF.classifications", colnames(x@meta.data), value = TRUE)
    colnames(x@meta.data)[match(df_col, colnames(x@meta.data))] <- "doublet_finder"
    x
  }, future.seed = TRUE)

  future::plan(sequential)

  # ---- Merge samples and remove predicted doublets ------------------------------
  combined_obj <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = accessions
  )
  combined_obj[["RNA"]] <- JoinLayers(combined_obj[["RNA"]])
  combined_obj[["RNA"]] <- split(combined_obj[["RNA"]], f = combined_obj$Seq_ID)
  DefaultAssay(combined_obj) <- "RNA"
  combined_obj <- subset(combined_obj, subset = doublet_finder == "Singlet")
  rm(seurat_list)

  message("Singlets retained: ", ncol(combined_obj))

  # ---- Normalization, dimensionality reduction, and integration ----------------
  combined_obj <- SCTransform(
    combined_obj,
    vars.to.regress = c("S.Score", "G2M.Score", "percent.mt")
  )
  combined_obj <- RunPCA(combined_obj)

  pct <- combined_obj[["pca"]]@stdev / sum(combined_obj[["pca"]]@stdev) * 100
  co1 <- which(cumsum(pct) > 90 & pct < 5)[1]
  co2 <- sort(which((pct[-length(pct)] - pct[-1]) > 0.1), decreasing = TRUE)[1] + 1
  num_pcs <- min(co1, co2)

  # Unintegrated representation retained for comparison/provenance.
  combined_obj <- FindNeighbors(combined_obj, dims = 1:num_pcs, reduction = "pca")
  combined_obj <- FindClusters(
    combined_obj,
    resolution = 0.1, cluster.name = "unintegrated_clusters"
  )
  combined_obj <- RunUMAP(
    combined_obj,
    dims = 1:num_pcs, reduction = "pca",
    reduction.name = "umap.unintegrated"
  )

  # Reproduce the three integration approaches evaluated in the original analysis.
  combined_obj <- IntegrateLayers(
    combined_obj,
    method = CCAIntegration, orig.reduction = "pca",
    new.reduction = "integrated.cca", normalization.method = "SCT"
  )
  combined_obj <- IntegrateLayers(
    combined_obj,
    method = RPCAIntegration, orig.reduction = "pca",
    new.reduction = "integrated.rpca", normalization.method = "SCT"
  )
  combined_obj <- IntegrateLayers(
    combined_obj,
    method = HarmonyIntegration, orig.reduction = "pca",
    new.reduction = "harmony", normalization.method = "SCT"
  )

  for (reduction in c("integrated.cca", "integrated.rpca", "harmony")) {
    umap_name <- switch(reduction,
      "integrated.cca" = "umap.cca",
      "integrated.rpca" = "umap.rpca",
      "harmony" = "umap.harmony"
    )
    cluster_name <- paste0(reduction, "_clusters")

    combined_obj <- RunUMAP(
      combined_obj,
      dims = 1:num_pcs, reduction = reduction,
      reduction.name = umap_name
    )
    combined_obj <- FindNeighbors(combined_obj, dims = 1:num_pcs, reduction = reduction)
    combined_obj <- FindClusters(
      combined_obj,
      resolution = 0.1, cluster.name = cluster_name
    )
  }

  # RPCA clustering was selected for downstream cell-type annotation.
  Idents(combined_obj) <- "integrated.rpca_clusters"

  cluster_labels <- c(
    "Urothelial cell", "T cell", "Endothelial", "MyoCAF", "iCAF", "Myeloid",
    "Proliferating cell", "B cell", "Differentiating plasma cell", "Mast cell",
    "Mature plasma cell"
  )
  names(cluster_labels) <- levels(combined_obj)
  combined_obj <- RenameIdents(combined_obj, cluster_labels)
  combined_obj$cell_type <- factor(
    as.character(Idents(combined_obj)),
    levels = cluster_labels
  )
  Idents(combined_obj) <- "cell_type"

  # ---- Add sample annotations from the Chen dataset ----------------------------
  meta <- combined_obj@meta.data |>
    tibble::rownames_to_column("cell_id") |>
    left_join(sample_metadata, by = "Seq_ID")

  if (anyNA(meta$Tumor_type))
    stop("Clinical metadata could not be matched to all SRA accessions.")

  meta$Tumor_type <- factor(meta$Tumor_type, levels = c("Normal", "NMIBC", "MIBC"))
  meta$Tumor_grade <- factor(
    meta$Tumor_grade,
    levels = c("Low grade", "High grade", "Adj. low grade", "Adj. high grade")
  )
  meta$Tumor_stage <- factor(meta$Tumor_stage, levels = c("NMIBC", "MIBC"))

  rownames(meta) <- meta$cell_id
  meta$cell_id <- NULL
  combined_obj@meta.data <- meta[colnames(combined_obj), , drop = FALSE]

  # ---- Save --------------------------------------------------------------------
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  saveRDS(combined_obj, out_file)
  message("Saved annotated Seurat object: ", out_file)
  invisible(combined_obj)
}
