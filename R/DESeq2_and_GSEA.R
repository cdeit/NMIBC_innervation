# ==============================================================================
# DESeq2 differential expression and gene set enrichment analysis
# ==============================================================================

library(DESeq2)
library(dplyr)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(stringr)
library(ggbreak)
library(ggplotify)
library(patchwork)

# ==============================================================================
# HELPERS
# ==============================================================================

# Null-coalescing helper for optional config values.
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Assign STaN-high and STaN-low groups by quantile
#'
#' Applies the configured quantile threshold separately within each cohort.
#' Samples between the low and high quantiles are excluded.
.apply_quantile_threshold <- function(pheno, cfg, cond = cfg$deseq2_condition) {
  q <- cfg$quantile_threshold
  if (is.null(q) || !grepl("^nerve_group", cond)) return(pheno)

  cont_col <- sub("^nerve_group", "nerve_obj_area_percent", cond)
  if (!cont_col %in% colnames(pheno)) return(pheno)

  cohorts <- if ("Cohort" %in% colnames(pheno)) unique(pheno$Cohort) else "All"

  out <- lapply(cohorts, function(cohort) {
    x <- if ("Cohort" %in% colnames(pheno)) pheno[pheno$Cohort == cohort, , drop = FALSE] else pheno
    qs <- quantile(x[[cont_col]], probs = c(q, 1 - q), na.rm = TRUE)

    group <- rep(NA_character_, nrow(x))
    group[x[[cont_col]] <= qs[1]] <- "Low"
    group[x[[cont_col]] >= qs[2]] <- "High"
    x[[cond]] <- factor(group, levels = c("Low", "High"))
    x[!is.na(x[[cond]]), , drop = FALSE]
  })

  do.call(rbind, out)
}


#' Format gene-set descriptions for plotting
.clean_gsea_descriptions <- function(x, gene_anns = NULL) {
  keep_lower <- c(
    "a", "the", "of", "and", "in", "on", "for", "with",
    "at", "by", "to", "from", "as", "is", "are"
  )
  lower_regex <- paste0("\\b(", paste(keep_lower, collapse = "|"), ")\\b")
  custom_cap <- c(
    "Mtorc1" = "mTORC1", "Jak" = "JAK", "Dna" = "DNA", "Mrna" = "mRNA",
    "Mirna" = "miRNA", "Rna" = "RNA", "Rnas" = "RNAs", "Trna" = "tRNA",
    "Atp" = "ATP", "Dntp" = "dNTP"
  )

  x <- x |>
    gsub("_", " ", x = _) |>
    stringr::str_to_title() |>
    stringr::str_replace_all(stringr::regex(lower_regex, ignore_case = TRUE), tolower) |>
    stringr::str_replace_all(custom_cap) |>
    stringr::str_replace("^([a-z])", toupper) |>
    stringr::str_replace_all("\\b([b-z])\\b", toupper)

  # Restore all-caps gene symbols when annotations are supplied.
  if (!is.null(gene_anns) && all(c("GENE_NAME", "GENE_TYPE") %in% colnames(gene_anns))) {
    genes <- gene_anns |>
      dplyr::filter(grepl("protein_coding", GENE_TYPE)) |>
      dplyr::pull(GENE_NAME)

    if (length(genes)) {
      gene_pattern <- paste0("\\b(", paste(genes, collapse = "|"), ")\\b")
      x <- stringr::str_replace_all(
        x, stringr::regex(gene_pattern, ignore_case = TRUE), toupper
      )
    }
  }

  x
}


# ==============================================================================
# DESeq2
# ==============================================================================

#' Run DESeq2 on combined RNA-seq cohorts
#'
#' Combines Cohorts A and B, applies STaN quantile stratification, filters genes
#' by configured biotype/expression thresholds, and fits a covariate-adjusted
#' DESeq2 model. The reported contrast is STaN-high vs STaN-low.
run_combined_deseq2 <- function(pheno_data, a_counts, b_counts, gene_anns, cfg, save = TRUE) {
  message("Running DESeq2 on Cohorts A & B combined...")
  cond <- cfg$deseq2_condition

  # Output path.
  out_dir <- file.path(cfg$path, cfg$paths$results$bulk_rnaseq, "tables")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_csv <- file.path(
    out_dir,
    paste0(
      "DESeq2_Results_Cohorts_AB_adj_",
      paste(cfg$deseq2_covariates, collapse = "_"), ".csv"
    )
  )

  # Stratify samples and remove rows missing requested covariates.
  pheno_data <- pheno_data[!is.na(pheno_data[["RNA-seq"]]), , drop = FALSE]
  if (!is.numeric(pheno_data[[cond]])) {
    pheno <- .apply_quantile_threshold(pheno_data, cfg, cond)
  } else (pheno <- pheno_data[[cond]])

  for (cov in cfg$deseq2_covariates)
    pheno <- pheno[!is.na(pheno[[cov]]), , drop = FALSE]

  pheno[[cond]] <- relevel(factor(pheno[[cond]]), ref = "Low")


  # Combine count matrices using genes shared across cohorts.
  if (cfg$deseq2_cohorts == "A") {
    b_counts <- NULL
  } else if (cfg$deseq2_cohorts == "B") {
    a_counts <- NULL
  } else if (cfg$deseq2_cohorts != "AB") {
    stop("Invalid value for cfg$deseq2_cohorts. Must be 'A', 'B', or 'AB'.")
  }

  if (is.null(a_counts)) {
    counts <- b_counts
  } else if (is.null(b_counts)) {
    counts <- a_counts
  } else {
    common_genes <- intersect(rownames(a_counts), rownames(b_counts))
    counts <- cbind(
      a_counts[common_genes, , drop = FALSE],
      b_counts[common_genes, , drop = FALSE]
    )
  }

  if (is.null(counts)) stop("Both a_counts and b_counts are NULL.")

  # Retain requested gene biotypes.
  valid_genes <- gene_anns$GENE_NAME[gene_anns$GENE_TYPE %in% cfg$gene_type_filter]
  counts <- counts[rownames(counts) %in% valid_genes, , drop = FALSE]

  # Align phenotype rows and count-matrix columns.
  pts <- intersect(pheno$`RNA-seq`, colnames(counts))
  pheno <- pheno[match(pts, pheno$`RNA-seq`), , drop = FALSE]
  print(table(pheno[[cond]]))

  for (cov in cfg$deseq2_covariates) {
    pheno <- pheno[!is.na(pheno[[cov]]) & pheno[[cov]] != "Missing" &
      pheno[[cov]] != "" & pheno[[cov]] != "Unknown", , drop = FALSE]
  }

  counts <- counts[, pts, drop = FALSE]

  # Remove low-expression genes.
  keep <- rowSums(counts >= cfg$min_count) >= cfg$min_samples
  counts <- counts[keep, , drop = FALSE]

  # Fit covariate-adjusted DESeq2 model.
  cov_str <- paste(cfg$deseq2_covariates, collapse = " + ")
  design <- as.formula(if (nzchar(cov_str))
    paste("~", cov_str, "+", cond) else paste("~", cond))

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = pheno,
    design = design
  )
  dds <- DESeq2::DESeq(dds)

  # Differential expression using dynamic levels
  exp_level <- levels(pheno[[cond]])[2]
  ref_level <- cfg$deseq2_ref

  res <- as.data.frame(DESeq2::results(dds, contrast = c(cond, exp_level, ref_level)))
  res$Gene <- rownames(res)

  # Add supplied gene annotations and classify significant DEGs.
  ann <- gene_anns[!duplicated(gene_anns$GENE_NAME), , drop = FALSE]
  res <- res |>
    dplyr::left_join(ann, by = c("Gene" = "GENE_NAME")) |>
    dplyr::relocate(Gene, .before = baseMean) |>
    dplyr::arrange(padj) |>
    dplyr::mutate(
      Significance = dplyr::case_when(
        padj < cfg$deg_padj_cutoff &
          log2FoldChange > cfg$deg_log2fc_cutoff ~ "Upregulated",
        padj < cfg$deg_padj_cutoff &
          log2FoldChange < -cfg$deg_log2fc_cutoff ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    ) |>
    as.data.frame()

  write.csv(res, out_csv, row.names = FALSE)
  message("Saved DESeq2 results to ", out_csv)

  list(res = res, dds = dds, cond = cond, exp_level = exp_level, ref_level = ref_level)
}


#' Plot DESeq2 volcano plot
#'
#' Labels either the strongest significant genes, equal numbers in each
#' direction, or a user-specified gene set according to label_method.
plot_deg_volcano <- function(res, cfg, save = TRUE, print = TRUE,
                             label_n = 50, label_method = "top", label_custom_genes = NULL) {
  out_dir <- file.path(cfg$path, cfg$paths$results$bulk_rnaseq, "figures")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_pdf <- file.path(
    out_dir,
    paste0(
      "DESeq2_DEG_Volcano_Cohorts_AB_adj_",
      paste(cfg$deseq2_covariates, collapse = "_"), ".pdf"
    )
  )

  degs <- res
  if (!"Significance" %in% colnames(degs)) {
    degs <- degs |>
      dplyr::mutate(
        Significance = dplyr::case_when(
          padj < cfg$deg_padj_cutoff &
            log2FoldChange > cfg$deg_log2fc_cutoff ~ "Upregulated",
          padj < cfg$deg_padj_cutoff &
            log2FoldChange < -cfg$deg_log2fc_cutoff ~ "Downregulated",
          TRUE ~ "Not Significant"
        )
      )
  }
  degs$Significance <- factor(
    degs$Significance,
    levels = c("Upregulated", "Downregulated", "Not Significant")
  )

  p <- ggplot2::ggplot(
    degs,
    ggplot2::aes(log2FoldChange, -log10(padj), color = Significance)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::scale_color_manual(values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "gray70"
    )) +
    ggplot2::geom_vline(
      xintercept = c(-cfg$deg_log2fc_cutoff, cfg$deg_log2fc_cutoff),
      linetype = "dashed", alpha = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(cfg$deg_padj_cutoff),
      linetype = "dashed", alpha = 0.5
    ) +
    ggplot2::labs(
      title = "STaN-High vs STaN-Low",
      x = "log2 Fold Change", y = "-log10 FDR",
      caption = sprintf(
        "Total Evaluated Genes: %d\nSignificantly Upregulated: %d\nSignificantly Downregulated: %d",
        nrow(degs),
        sum(degs$Significance == "Upregulated", na.rm = TRUE),
        sum(degs$Significance == "Downregulated", na.rm = TRUE)
      )
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.background = ggplot2::element_rect(fill = "white", color = NA) # ,
      # aspect.ratio = 1
    )

  # Select genes to label.
  top_df <- NULL
  if (!is.null(label_n) && label_n > 0) {
    method <- label_method %||% "top"

    if (method == "symmetric") {
      n_each <- max(1, ceiling(label_n / 2))
      top_df <- dplyr::bind_rows(
        degs |>
          dplyr::filter(
            !is.na(stat), log2FoldChange > cfg$deg_log2fc_cutoff,
            padj < cfg$deg_padj_cutoff
          ) |>
          dplyr::arrange(dplyr::desc(stat)) |>
          utils::head(n_each),
        degs |>
          dplyr::filter(
            !is.na(stat), log2FoldChange < -cfg$deg_log2fc_cutoff,
            padj < cfg$deg_padj_cutoff
          ) |>
          dplyr::arrange(stat) |>
          utils::head(n_each)
      )
    } else if (method == "custom") {
      top_df <- degs |>
        dplyr::filter(Gene %in% label_custom_genes)
    } else {
      top_df <- degs |>
        dplyr::filter(
          !is.na(stat), padj < cfg$deg_padj_cutoff,
          abs(log2FoldChange) > cfg$deg_log2fc_cutoff
        ) |>
        dplyr::arrange(dplyr::desc(abs(stat))) |>
        utils::head(label_n)
    }

    if (nrow(top_df)) {
      p <- p + ggrepel::geom_text_repel(
        data = top_df,
        ggplot2::aes(label = Gene),
        size = 3, color = "black", fontface = "bold",
        min.segment.length = 0.25, max.overlaps = 15
      )
    }
  }

  if (save) {
    ggplot2::ggsave(out_pdf, p, width = 7, height = 7, limitsize = FALSE)
    message("Saved volcano plot to ", out_pdf)
  }

  if (print) {
    p
  }
}


# ==============================================================================
# GSEA
# ==============================================================================

#' Collapse redundant gene sets by leading-edge Jaccard similarity
.collapse_gsea_jaccard <- function(cp_obj, cutoff) {
  df <- cp_obj@result
  if (nrow(df) <= 1) return(cp_obj)

  df <- df[order(df$p.adjust), , drop = FALSE]
  genes <- lapply(df$core_enrichment, \(x) strsplit(as.character(x), "/")[[1]])

  keep <- 1L
  for (i in 2:length(genes)) {
    genes_i <- genes[[i]]
    is_redundant <- FALSE

    for (j in keep) {
      genes_j <- genes[[j]]
      intersect_len <- length(intersect(genes_i, genes_j))
      union_len <- length(union(genes_i, genes_j))

      if (union_len > 0 && (intersect_len / union_len) > cutoff) {
        is_redundant <- TRUE
        break
      }
    }

    if (!is_redundant) keep <- c(keep, i)
  }

  cp_obj@result <- df[keep, , drop = FALSE]
  cp_obj
}


#' Optionally collapse semantically redundant GO Biological Process terms
#'
#' Semantic simplification is applied only when the GMT term IDs are GO accessions.
#' Other GMTs, including GOBP files containing only text pathway names, are left
#' unchanged and still receive the universal Jaccard redundancy filter.
.simplify_gobp <- function(cp_obj, dataset_name, cfg) {
  is_gobp <- grepl("GOBP|GO[_.-]?BP|GO.*BIOLOGICAL.*PROCESS",
    dataset_name,
    ignore.case = TRUE
  )
  if (!is_gobp || !nrow(cp_obj@result)) return(cp_obj)

  go_ids <- cp_obj@result$ID
  if (!all(grepl("^GO:[0-9]{7}$", go_ids))) {
    message(
      "GO semantic collapse skipped for ", dataset_name,
      ": GMT terms are not GO accession IDs."
    )
    return(cp_obj)
  }

  sem_data <- GOSemSim::godata(
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = "BP",
    computeIC = FALSE
  )

  clusterProfiler::simplify(
    cp_obj,
    cutoff = cfg$semantic_cutoff,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = sem_data
  )
}


#' Run pre-ranked GSEA for one GMT file
#'
#' Uses the DESeq2 Wald statistic as the ranking metric, optionally performs GO
#' semantic simplification for GOBP GMTs, then removes remaining redundant terms
#' by leading-edge Jaccard similarity.

run_gsea_gmt <- function(res, gmt_file, dataset_name = NULL, cfg) {
  if (is.null(dataset_name))
    dataset_name <- tools::file_path_sans_ext(basename(gmt_file))

  gmt_path <- if (file.exists(gmt_file)) gmt_file else
    file.path(cfg$path, "data", "genesets", gmt_file)
  if (!file.exists(gmt_path)) stop("GMT file not found: ", gmt_path)

  out_dir <- file.path(cfg$path, cfg$paths$results$bulk_rnaseq, "tables")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_csv <- file.path(out_dir, paste0("GSEA_Results_", dataset_name, ".csv"))
  message("Running GSEA for ", dataset_name, "...")

  ranks <- res[!is.na(res$stat) & !is.na(res$Gene), c("Gene", "stat")]
  ranks <- ranks[!duplicated(ranks$Gene), , drop = FALSE]
  gene_ranks <- sort(stats::setNames(ranks$stat, ranks$Gene), decreasing = TRUE)

  term2gene <- clusterProfiler::read.gmt(gmt_path)

  set.seed(42)
  cp_obj <- clusterProfiler::GSEA(
    geneList = gene_ranks,
    TERM2GENE = term2gene,
    minGSSize = cfg$gsea_min_size,
    maxGSSize = cfg$gsea_max_size,
    pvalueCutoff = 1,
    eps = cfg$gsea_eps
  )

  # Guard against NULL if no gene sets pass thresholds
  if (is.null(cp_obj) || nrow(cp_obj@result) == 0) {
    message("No enriched gene sets found for ", dataset_name)
    return(NULL)
  }

  cp_obj <- .simplify_gobp(cp_obj, dataset_name, cfg)
  cp_obj <- .collapse_gsea_jaccard(cp_obj, cfg$jaccard_cutoff)

  cp_obj@result <- cp_obj@result |>
    dplyr::mutate(
      Significance = dplyr::if_else(
        p.adjust < cfg$gsea_padj_cutoff, "Significant", "Not Significant"
      ),
      Enrichment = factor(
        dplyr::if_else(NES > 0, "STaN-high", "STaN-low"),
        levels = c("STaN-low", "STaN-high")
      )
    )

  write.csv(cp_obj@result, out_csv, row.names = FALSE)
  message("Saved GSEA results to ", out_csv)
  cp_obj
}

#' Run all GMT gene-set analyses and associated plots
#'
#' Uses the named GMT collections in cfg$gmt_files and skips plots for empty results.
#' Returns GSEA and available plot objects in a named list.
run_all_gsea <- function(res, cfg, gene_anns = NULL,
                         geneset_dir = file.path(cfg$path, "data", "genesets"),
                         make_plots = TRUE, print = TRUE) {
  if (is.null(cfg$gmt_files) || length(cfg$gmt_files) == 0) {
    stop("No GMT files specified in cfg$gmt_files")
  }

  dataset_names <- names(cfg$gmt_files)
  gmt_files <- unname(unlist(cfg$gmt_files))
  message("Found ", length(gmt_files), " GMT file(s) in config: ",
          paste(dataset_names, collapse = ", "))
  out <- stats::setNames(vector("list", length(gmt_files)), dataset_names)

  for (i in seq_along(gmt_files)) {
    nm <- dataset_names[i]
    gmt_file <- file.path(geneset_dir, gmt_files[i])
    cp_obj <- run_gsea_gmt(res, gmt_file, nm, cfg)
    if (is.null(cp_obj)) {
      out[[nm]] <- list(gsea = NULL)
      next
    }
    out[[nm]] <- list(gsea = cp_obj)
    if (isTRUE(make_plots)) {
      out[[nm]]$nes_plot <- tryCatch(
        plot_nes_dotplot(cp_obj, nm, cfg, gene_anns = gene_anns),
        error = function(e) {
          message("NES dotplot skipped for ", nm, ": ", conditionMessage(e))
          NULL
        }
      )
      out[[nm]]$jaccard_plot <- tryCatch(
        plot_jaccard_heatmap(cp_obj, nm, cfg),
        error = function(e) {
          message("Jaccard heatmap skipped for ", nm, ": ", conditionMessage(e))
          NULL
        }
      )
    }
  }
  if (isTRUE(make_plots) && isTRUE(print)) {
    for (nm in dataset_names) {
      if (!is.null(out[[nm]]$nes_plot)) print(out[[nm]]$nes_plot)
      if (!is.null(out[[nm]]$jaccard_plot)) print(out[[nm]]$jaccard_plot)
    }
  }
  invisible(out)
}


# ==============================================================================
# GSEA PLOTS
# ==============================================================================

#' Plot normalized enrichment scores
#'
#' Shows the top positively and negatively enriched gene sets, with point color
#' representing FDR and point size representing leading-edge gene count.
plot_nes_dotplot <- function(cp_obj, dataset_name, cfg, save = TRUE, print = TRUE, gene_anns = NULL, label_width = 35) {
  out_dir <- file.path(cfg$path, cfg$paths$results$bulk_rnaseq, "figures")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (dataset_name %in% names(cfg$gmt_files)) {
    dataset_name_short <- dataset_name
  } else if (!is.null(cfg$gmt_files) && dataset_name %in% unlist(cfg$gmt_files)) {
    dataset_name_short <- names(
      cfg$gmt_files
    )[sapply(cfg$gmt_files, function(x) x == dataset_name)]
  } else {
    dataset_name_short <- gsub("_|\\.gmt\\.txt|\\.gmt$", " ", dataset_name)
  }

  out_pdf <- file.path(out_dir, paste0(
    "NES_Dotplot_simcutoff", cfg$jaccard_cutoff,
    "_top", cfg$top_gsets_plot_n, "_", dataset_name, ".pdf"
  ))

  is_cell_sig <- grepl("Panglao|cell", dataset_name, ignore.case = TRUE)
  title <- paste(
    dataset_name_short,
    if (is_cell_sig) "Cell Type Signature Enrichment" else "Gene Set Enrichment"
  )

  df <- cp_obj@result
  if (!nrow(df)) stop("No results to plot for ", dataset_name)

  df$Description <- .clean_gsea_descriptions(df$Description, gene_anns)

  # Optionally limit plotting to gene sets meeting both FDR and NES thresholds.
  if (isTRUE(cfg$sig_gsets_only)) {
    df <- df[
      df$p.adjust < cfg$gsea_padj_cutoff &
        abs(df$NES) >= cfg$gsea_nes_cutoff, ,
      drop = FALSE
    ]
  }
  if (!nrow(df)) stop("No gene sets pass the plotting thresholds for ", dataset_name)

  # Select top positive and negative gene sets by FDR.
  up <- df[df$NES > cfg$gsea_nes_cutoff, , drop = FALSE]
  dn <- df[df$NES < -cfg$gsea_nes_cutoff, , drop = FALSE]
  up <- utils::head(up[order(up$p.adjust), , drop = FALSE], cfg$top_gsets_plot_n)
  dn <- utils::head(dn[order(dn$p.adjust), , drop = FALSE], cfg$top_gsets_plot_n)
  df_plot <- rbind(up, dn)

  if (!nrow(df_plot)) stop("No gene sets exceed the NES cutoff for ", dataset_name)

  df_plot$Description <- factor(
    df_plot$Description,
    levels = df_plot$Description[order(df_plot$NES)]
  )
  df_plot$setSize <- vapply(
    strsplit(as.character(df_plot$core_enrichment), "/"),
    length, integer(1)
  )

  # Dynamic dimensions preserve approximately constant row/panel sizing.
  n_wrapped <- sum(nchar(as.character(df_plot$Description)) > label_width)
  plot_h <- max(3, 1.5 + nrow(df_plot) * 0.22 + n_wrapped * 0.66)
  plot_w <- 4.25 + label_width * 0.06

  caption_hjust <- 0.9262795 # Preserves original caption alignment.
  leg_h <- 6 * 0.22
  leg_w <- 0.85
  leg_x <- 1 - 0.3 / plot_w
  leg_y <- 0.7 / plot_h

  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(NES, Description, color = p.adjust, size = setSize)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "red", high = "blue") +
    ggplot2::labs(y = NULL, color = "FDR", size = "Gene Count", title = title) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = 0.06),
      breaks = unique(c(
        round(min(df_plot$NES, na.rm = TRUE) - 0.5), -2, -1, 1, 2,
        round(max(df_plot$NES, na.rm = TRUE) + 0.5)
      )),
      limits = c(
        round(min(df_plot$NES, na.rm = TRUE) - 0.5),
        round(max(df_plot$NES, na.rm = TRUE) + 0.5)
      )
    ) +
    ggplot2::theme(
      text = element_text(family = "Helvetica"),
    ) +
    ggplot2::scale_y_discrete(
      labels = \(x) ifelse(
        nchar(x) > label_width,
        stringr::str_wrap(x, width = label_width),
        x
      )
    )

  # Build legend separately so its dimensions remain stable across plots.
  size_breaks <- c(25, 50, 100, 200, 300)
  size_breaks <- size_breaks[size_breaks <= max(df_plot$setSize)]
  if (!length(size_breaks)) size_breaks <- pretty(df_plot$setSize, n = 3)

  p_leg <- p +
    ggplot2::scale_size_continuous(
      breaks = utils::tail(size_breaks, if (nrow(df_plot) <= 10) 2 else 5)
    ) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(
        size = 11, margin = ggplot2::margin(b = 1, t = 6.5, unit = "pt")
      ),
      legend.position = "inside",
      legend.location = "full",
      legend.key = ggplot2::element_rect(fill = NA, color = NA),
      legend.text = ggplot2::element_text(size = 11),
      legend.key.size = grid::unit(0.45, "cm"),
      legend.spacing.y = grid::unit(0.15, "cm"),
      legend.margin = ggplot2::margin(t = -0.2, unit = "cm"),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 0, l = 10, unit = "pt")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(barheight = 3.5, order = 1),
      size = ggplot2::guide_legend(order = 2)
    )
  leg <- ggpubr::get_legend(p_leg)

  # Remove the uninformative central NES interval to emphasize enriched terms.
  p <- p +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "plain", hjust = 0.75),
      plot.title.position = "plot",
      axis.title = ggplot2::element_blank(),
      axis.line.x.top = ggplot2::element_blank(),
      axis.ticks.x.top = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_blank(),
      axis.text.y.left = ggplot2::element_text(size = 13),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 7.5, l = 10, unit = "pt")
    ) +
    ggbreak::scale_x_break(c(-0.75, 0.75), space = 0.3, scales = "fixed", expand = TRUE)

  p <- ggplotify::as.ggplot(print(p))
  p_final <- patchwork::wrap_elements(p) +
    patchwork::plot_annotation(
      caption = "Normalized Enrichment Score",
      theme = ggplot2::theme(
        text = element_text(family = "Helvetica"),
        plot.caption.position = "plot",
        plot.caption = ggplot2::element_text(
          size = 12,
          margin = ggplot2::margin(t = -0.45, unit = "cm"),
          hjust = caption_hjust
        )
      )
    )

  patch <- cowplot::ggdraw(ggplotify::as.ggplot(p_final)) +
    cowplot::draw_grob(
      leg,
      x = leg_x, y = leg_y,
      height = leg_h / plot_h,
      width = leg_w / plot_w,
      hjust = 1, vjust = 0
    )

  if (isTRUE(save)) {
    cowplot::ggsave2(
      patch,
      filename = out_pdf,
      width = plot_w, height = plot_h, units = "in"
    )

    message("Saved NES dotplot to ", out_pdf)
  }

  message("NES < 0 = STaN-low; NES > 0 = STaN-high")

  if (print) patch else invisible(patch)
}


#' Plot pairwise Jaccard overlap of leading-edge genes
#'
#' Displays the leading-edge gene overlap among the strongest positively and
#' negatively enriched pathways.
plot_jaccard_heatmap <- function(cp_obj, dataset_name, cfg, make_plot = TRUE, print = FALSE) {
  out_dir <- file.path(cfg$path, cfg$paths$results$bulk_rnaseq, "figures")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (dataset_name %in% names(cfg$gmt_files)) {
    dataset_name_short <- dataset_name
  } else if (!is.null(cfg$gmt_files) && dataset_name %in% unlist(cfg$gmt_files)) {
    dataset_name_short <- names(
      cfg$gmt_files
    )[sapply(cfg$gmt_files, function(x) x == dataset_name)]
  } else {
    dataset_name_short <- gsub("_|\\.gmt\\.txt|\\.gmt$", " ", dataset_name)
  }

  out_pdf <- file.path(out_dir, paste0(
    "Jaccard_Heatmap_simcutoff", cfg$jaccard_cutoff,
    "_top", cfg$top_gsets_plot_n, "_", dataset_name_short, ".pdf"
  ))

  df <- cp_obj@result
  up <- df[df$NES > 0, , drop = FALSE]
  dn <- df[df$NES < 0, , drop = FALSE]

  up_pw <- utils::head(
    up$Description[order(up$NES, decreasing = TRUE)],
    cfg$top_gsets_plot_n
  )
  dn_pw <- utils::head(
    dn$Description[order(dn$NES)],
    cfg$top_gsets_plot_n
  )
  pathways <- c(up_pw, dn_pw)

  if (length(pathways) < 2)
    stop("Insufficient pathways for Jaccard overlap matrix in ", dataset_name)

  edges <- lapply(pathways, function(pw) {
    x <- df$core_enrichment[df$Description == pw][1]
    unlist(strsplit(as.character(x), "/"))
  })

  n <- length(pathways)
  j_mat <- matrix(0, n, n, dimnames = list(pathways, pathways))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      union_n <- length(union(edges[[i]], edges[[j]]))
      j_mat[i, j] <- if (i == j) 1 else if (!union_n) 0 else
        length(intersect(edges[[i]], edges[[j]])) / union_n
    }
  }

  melted <- expand.grid(
    Pathway1 = pathways,
    Pathway2 = rev(pathways),
    stringsAsFactors = FALSE
  )
  melted$Jaccard <- mapply(
    \(x, y) j_mat[x, y],
    melted$Pathway1, melted$Pathway2
  )
  melted$Pathway1 <- factor(melted$Pathway1, levels = pathways)
  melted$Pathway2 <- factor(melted$Pathway2, levels = rev(pathways))

  p <- ggplot2::ggplot(
    melted,
    ggplot2::aes(Pathway1, Pathway2, fill = Jaccard)
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0.5, name = "Jaccard Similarity"
    ) +
    ggplot2::labs(
      title = paste("Jaccard Similarity of Leading Edge Genes -", dataset_name),
      subtitle = "Red labels = STaN-high; Blue labels = STaN-low"
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, size = 10,
        color = ifelse(levels(melted$Pathway1) %in% up_pw, "darkred", "darkblue")
      ),
      axis.text.y = ggplot2::element_text(
        size = 10,
        color = ifelse(levels(melted$Pathway2) %in% up_pw, "darkred", "darkblue")
      ),
      axis.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "plain", size = 20),
      aspect.ratio = 1
    )

  max_char <- max(nchar(pathways))
  plot_w <- 4 + n * 0.2 + max_char * 0.06
  plot_h <- 3 + n * 0.2 + max_char * 0.06

  if (make_plot) {
    ggplot2::ggsave(
      out_pdf, p,
      width = plot_w, height = plot_h, limitsize = FALSE
    )
    message("Saved Jaccard heatmap to ", out_pdf)
  }

  if (print) p else invisible(p)
}
