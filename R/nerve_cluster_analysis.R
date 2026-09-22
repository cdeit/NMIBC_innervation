#' ============================================================
# --- Nerve (STaN) clusters -----------------------------------
#' ============================================================


##  --- Helper functions --------------------------------------

.split_im3s <- function(x) {
  # Split comma-separated image IDs and return unique values
  unique(trimws(unlist(strsplit(stats::na.omit(x), ",\\s*"))))
}

.clean_marker_names <- function(d) {
  # Standardize marker names exported from inForm
  names(d) <- sub("^Object ", "", names(d))
  names(d) <- sub(" \\(.*$| Mean.*$", "", names(d))
  names(d)[names(d) == "Autofluorescence"] <- "AF"
  d
}

.cap_minmax <- function(x, q) {
  # Cap extreme values at the specified quantile, then rescale to 0-1
  x <- pmin(x, stats::quantile(x, q, na.rm = TRUE))
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

.save_subcluster_plot <- function(p, filename, cfg, width, height) {
  # Build output path and save plot
  out <- file.path(cfg$path, cfg$paths$results$nerve_clusters, "figures", filename)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out, p, width = width, height = height, units = "in")
  message("Saved: ", normalizePath(out, mustWork = FALSE))
}

.prep_subcluster_patients <- function(data, cfg, composition) {
  s <- cfg$nerve_subclusters

  # Restrict to samples with sufficient stromal area
  out <- dplyr::filter(
    data,
    .data[[s$stroma_area_col]] > s$min_stroma_area
  )

  # For composition plots, require appreciable total STaN abundance
  if (composition)
    out <- dplyr::filter(
      out,
      .data[[s$total_nerve_col]] > s$composition_min_total
    )

  # Standardize historical abundance-column naming
  names(out) <- sub(
    "final_obj_area_percent", "nerve_obj_area_percent",
    names(out), fixed = TRUE
  )

  out
}

.subcluster_long <- function(df, outcome_col, s) {
  # Convert subcluster abundance columns to long format
  df %>%
    dplyr::select(
      dplyr::all_of(c(outcome_col, unname(s$subcluster_cols)))
    ) %>%
    tidyr::pivot_longer(
      -dplyr::all_of(outcome_col),
      names_to = "raw",
      values_to = "NerveAreaPct"
    ) %>%

    # Map raw abundance columns to configured subtype names
    dplyr::mutate(
      Outcome = factor(.data[[outcome_col]]),
      Subtype = factor(
        names(s$subcluster_cols)[match(raw, s$subcluster_cols)],
        levels = s$subcluster_levels
      )
    ) %>%
    dplyr::filter(
      !is.na(Outcome),
      !is.na(NerveAreaPct),
      !is.na(Subtype)
    )
}


##  --- Main functions --------------------------------------


#' Plot marker enrichment by nerve subcluster
#'
#' Filters the stored cluster-versus-global marker deltas to the configured
#' markers, plots one facet per nerve subcluster, and exports a PDF.
#'
#' @param res_nerve_subclust Nerve-subclustering result containing
#'   `top_markers_per_nerve_cluster_list[[cfg$nerve_subclusters$marker_result]]`.
#'   That table must contain the configured marker, delta, and cluster columns.
#' @param cfg Analysis configuration containing `path`, `nerve_subclusters`, and
#'   `color_pals$nerve_subclusters`; see `config.R`.
#' @param filename PDF filename written beneath the configured output directory.
#'
#' @return Invisibly returns the `ggplot` object. The PDF is the file output.
#'
#' @examples
#' p_enrich <- plot_subcluster_marker_enrichment(res_nerve_subclust, cfg)
plot_subcluster_marker_enrichment <- function(res_nerve_subclust, cfg,
                                              filename = "subcluster_marker_enrichment.pdf") {

  s <- cfg$nerve_subclusters; d <- res_nerve_subclust$top_markers_per_nerve_cluster_list[[s$marker_result]]
  req <- c(s$marker_col, s$delta_col, s$cluster_col); if (!all(req %in% names(d))) stop("Missing marker-enrichment columns: ", paste(setdiff(req, names(d)), collapse = ", "))
  d <- dplyr::filter(d, .data[[s$marker_col]] %in% s$markers_keep) %>%
    dplyr::mutate(!!s$marker_col := factor(.data[[s$marker_col]], levels = s$markers_keep))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[s$marker_col]], y = .data[[s$delta_col]], fill = .data[[s$cluster_col]])) +
    ggplot2::geom_col() + ggplot2::geom_hline(yintercept = 0, color = "grey20") +
    ggplot2::facet_wrap(stats::as.formula(paste("~", s$cluster_col)), scales = "free_y", ncol = 1) +
    ggplot2::scale_fill_manual(values = cfg$color_pals$nerve_subclusters) +
    ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    ggplot2::labs(title = "Enriched markers per cluster", x = NULL, y = "Expression Delta (Cluster - Global)") +
    ggplot2::coord_cartesian(clip = "off") + cowplot::theme_cowplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1.1), panel.spacing.y = grid::unit(0.08, "cm")) +
    ggplot2::guides(fill = "none")
  .save_subcluster_plot(p, filename, cfg, width = 5, height = 8.5); invisible(p)
}

#' Generate and plot the nerve-object UMAP
#'
#' Generates a UMAP from the originally scaled object-level marker matrix,
#' restricts the displayed objects to the final image cohort, and plots both
#' subcluster identity and scaled marker MPI. Marker MPI used for coloring is
#' capped and min-max scaled independently; the UMAP input retains the original
#' scaling in `res$scaled_objs4pca_k`.
#'
#' @param res_nerve_subclust Nerve-subclustering result. `$nerve_objs` must
#'   contain the configured sample, object ID, cluster, and marker-mean columns.
#' @param res Original clustering result containing `$scaled_objs4pca_k`, with
#'   sample and object IDs plus every marker in `umap_input_markers`.
#' @param patient_data Final patient-level dataset containing `consolidated_im3s`
#'   (or the column configured by `im3_col`).
#' @param original_patient_data Patient-level dataset used when the object marker
#'   matrix was originally scaled; it must contain the same image-list column.
#' @param cfg Analysis configuration; see `config.R`.
#' @param filename PDF filename written beneath the configured output directory.
#'
#' @return Invisibly returns a list with `data`, `cluster_plot`,
#'   `expression_plots`, and `combined_plot`. The combined plot is saved as PDF.
#'
#' @examples
#' res_umap <- plot_nerve_umaps(res_nerve_subclust, res,
#'   patient_nerve_quant_07302025,
#'   patient_summary_nerve_quant_07302025_FINAL, cfg)
plot_nerve_umaps <- function(res_nerve_subclust, res, patient_data, original_patient_data, cfg,
                             filename = "nerve_object_umaps.pdf") {
  s <- cfg$nerve_subclusters; set.seed(s$umap_seed)
  final_im3s <- .split_im3s(patient_data[[s$im3_col]]); og_im3s <- .split_im3s(original_patient_data[[s$im3_col]])
  objs <- res_nerve_subclust$nerve_objs; marker_df <- objs %>%
    dplyr::filter(.data[[s$sample_col]] %in% final_im3s) %>%
    dplyr::select(dplyr::all_of(c(s$sample_col, s$object_col, s$cluster_col)), dplyr::contains("Mean")) %>%
    .clean_marker_names() %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(s$umap_plot_markers), .cap_minmax, q = s$marker_cap_quantile))

  d <- res$scaled_objs4pca_k %>%
    dplyr::filter(.data[[s$object_col]] %in% objs[[s$object_col]], .data[[s$sample_col]] %in% og_im3s) %>%
    .clean_marker_names() %>%
    dplyr::left_join(dplyr::select(objs, dplyr::all_of(c(s$object_col, s$cluster_col))), by = s$object_col)
  missing <- setdiff(s$umap_input_markers, names(d)); if (length(missing)) stop("Missing UMAP markers: ", paste(missing, collapse = ", "))
  emb <- uwot::umap(dplyr::select(d, dplyr::all_of(s$umap_input_markers)), n_neighbors = s$umap_n_neighbors,
                    min_dist = s$umap_min_dist, repulsion_strength = s$umap_repulsion, metric = "euclidean", verbose = TRUE)
  umap_df <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2], dplyr::select(d, dplyr::all_of(c(s$object_col, s$sample_col, s$cluster_col)))) %>%
    dplyr::filter(.data[[s$sample_col]] %in% final_im3s) %>%
    dplyr::left_join(dplyr::select(marker_df, dplyr::all_of(c(s$object_col, s$umap_plot_markers))), by = s$object_col)

  u <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color = .data[[s$cluster_col]])) +
    ggplot2::geom_point(size = 0.5, alpha = 0.9) + ggplot2::scale_color_manual(values = cfg$color_pals$nerve_subclusters) +
    ggplot2::coord_equal() + ggplot2::labs(title = "UMAP of Nerve Objects", color = "Subcluster") + cowplot::theme_cowplot() +
    ggplot2::theme(aspect.ratio = 1, axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
  expr <- lapply(s$umap_plot_markers, function(marker) ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color = .data[[marker]])) +
                   ggplot2::geom_point(size = 0.2, alpha = 0.8) +
                   ggplot2::scale_color_gradientn(name = "Scaled\nMPI", limits = c(0, 1), values = scales::rescale(c(0.15, 1)),
                                                  colours = cetcolor::cet_pal(9, name = "l8"), oob = scales::squish) +
                   ggplot2::coord_equal() + ggplot2::labs(title = marker) + cowplot::theme_cowplot() +
                   ggplot2::theme(aspect.ratio = 1, axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank()))

  patch <- patchwork::wrap_plots(u, patchwork::wrap_plots(expr, ncol = 2, guides = "collect"), ncol = 2) & ggplot2::theme(legend.position = "right")
  .save_subcluster_plot(patch, filename, cfg, width = 12, height = 8)
  invisible(list(data = umap_df, cluster_plot = u, expression_plots = expr, combined_plot = patch))
}

#' Plot mean subcluster abundance by a two-level outcome
#'
#' Filters the patient cohort using the configured stromal-area
#' criterion, plots mean subcluster abundance +/- SEM for a binary outcome, and
#' annotates per-subcluster adjusted P values.
#'
#' @param data Patient-level data containing the configured filtering columns,
#'   `outcome_col`, and all columns in `subcluster_cols`.
#' @param outcome_col Character scalar naming a two-level outcome, for example
#'   `"BCG_failure"` or `"Progression"`.
#' @param cfg Analysis configuration; see `config.R`.
#' @param filename PDF filename written beneath the configured output directory.
#'
#' @return Invisibly returns a list with `plot`, `statistics`, and long-form
#'   `data`. The plot is saved as PDF.
#'
#' @examples
#' res_bcg_abundance <- plot_cluster_abundance(
#'   patient_nerve_quant_07302025, "BCG_failure", cfg)
#' res_prog_abundance <- plot_cluster_abundance(
#'   patient_nerve_quant_07302025, "Progression", cfg)
plot_cluster_abundance <- function(data, cfg, outcome_col,
                                   filename = paste0("subcluster_abundance_by_", outcome_col, ".pdf"),
                                   bracket_offset = 0.01, bracket_tick = 0.002,
                                   outcome_colors = c("gray80", "gray30")) {
  s <- cfg$nerve_subclusters; df <- .prep_subcluster_patients(data, cfg, composition = FALSE)

  long <- .subcluster_long(df, outcome_col, s) %>%
    dplyr::mutate(Subtype = forcats::fct_reorder(Subtype, NerveAreaPct, median, na.rm = TRUE))
  if (nlevels(long$Outcome) != 2) stop("`", outcome_col, "` must have exactly two observed levels.")
  bars <- long %>% dplyr::group_by(Subtype, Outcome) %>%
    dplyr::summarise(mean_y = mean(NerveAreaPct), se_y = stats::sd(NerveAreaPct) / sqrt(dplyr::n()), .groups = "drop") %>%
    dplyr::group_by(Subtype) %>% dplyr::summarise(bar_top = max(mean_y + se_y, na.rm = TRUE), .groups = "drop")

  stats_df <- long %>% dplyr::group_by(Subtype) %>%
    {if (s$abundance_test == "t") rstatix::t_test(., NerveAreaPct ~ Outcome) else rstatix::wilcox_test(., NerveAreaPct ~ Outcome)} %>%
    dplyr::ungroup() %>% rstatix::adjust_pvalue(method = s$p_adjust) %>% dplyr::left_join(bars, by = "Subtype")
  yr <- diff(range(long$NerveAreaPct)); if (!is.finite(yr) || yr == 0) yr <- max(long$NerveAreaPct)
  stats_df <- stats_df %>% dplyr::mutate(p_label = dplyr::if_else(p.adj < 0.0001, "<0.0001", sprintf("%.3f", p.adj)), x = as.numeric(Subtype),
                                         xstart = x - 0.2, xend = x + 0.2, y = bar_top + bracket_offset * yr, tick = bracket_tick * yr)
  pd <- ggplot2::position_dodge(0.8)
  p <- ggplot2::ggplot(long, ggplot2::aes(Subtype, NerveAreaPct, fill = Outcome)) +
    ggplot2::stat_summary(fun = mean, geom = "col", position = pd, width = 0.7) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 18, size = 3, position = pd, show.legend = FALSE) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", position = pd, width = 0.2) +
    ggplot2::geom_segment(data = stats_df, ggplot2::aes(x = xstart, xend = xstart, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = 0.3) +
    ggplot2::geom_segment(data = stats_df, ggplot2::aes(x = xstart, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
    ggplot2::geom_segment(data = stats_df, ggplot2::aes(x = xend, xend = xend, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = 0.3) +
    ggplot2::geom_text(data = stats_df, ggplot2::aes(x, y + tick, label = p_label), inherit.aes = FALSE, size = 3, vjust = 0) +
    ggplot2::scale_fill_manual(values = outcome_colors) + ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = "Nerve subcluster", y = "Mean STaN+ area (%)", fill = gsub("_", " ", outcome_col)) + cowplot::theme_cowplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.margin = ggplot2::margin(5.5, 20, 10, 5.5)) +
    ggplot2::guides(fill = ggplot2::guide_legend(position = "top"))
  .save_subcluster_plot(p, filename, cfg, 6.5, 4.75); invisible(list(plot = p, statistics = stats_df, data = long))
}


#' Plot relative or absolute subcluster composition by outcome
#'
#' Summarizes each subcluster across patients within outcome groups and produces
#' either 100% stacked composition bars or absolute stacked abundance bars.
#' Patients must pass the stromal-area and minimum-total-STaN
#' filters configured in `cfg`.
#'
#' @param data Patient-level data containing the configured filtering columns,
#'   total nerve-abundance column, `outcome_col`, and all `subcluster_cols`.
#' @param outcome_col Character scalar naming the grouping outcome.
#' @param cfg Analysis configuration; see `config.R`.
#' @param relative If `TRUE`, normalize each outcome bar to 100%; otherwise show
#'   the absolute stacked summary.
#' @param stat Summary across patients: `"mean"` or `"median"`.
#' @param filename PDF filename written beneath the configured output directory.
#'
#' @return Invisibly returns a list with `plot`, summarized `summary`, and
#'   long-form `data`. The plot is saved as PDF.
#'
#' @examples
#' res_bcg_composition <- plot_subcluster_composition(
#'   patient_nerve_quant_07302025, "BCG_failure", cfg)
#' res_prog_absolute <- plot_cluster_composition(
#'   patient_nerve_quant_07302025, "Progression", cfg, relative = FALSE)
plot_cluster_composition <- function(data, cfg, outcome_col, relative = TRUE, stat = "mean",
                                     filename = paste0("subcluster_composition_by_", outcome_col, ".pdf")) {
  s <- cfg$nerve_subclusters; df <- .prep_subcluster_patients(data, cfg, composition = TRUE); long <- .subcluster_long(df, outcome_col, s)
  sum_df <- long %>% dplyr::group_by(Outcome, Subtype) %>%
    dplyr::summarise(value = if (stat == "mean") mean(NerveAreaPct) else median(NerveAreaPct), .groups = "drop") %>%
    # dplyr::mutate(Subtype = factor(Subtype, levels = s$subcluster_levels))

    # reorder subtype by overall contribution (optional but usually nicer)
    dplyr::group_by(Subtype) %>%
    dplyr::mutate(overall = sum(value, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Subtype = forcats::fct_reorder(Subtype, overall))

  p <- ggplot2::ggplot(sum_df, ggplot2::aes(Outcome, value, fill = Subtype)) +
    ggplot2::geom_col(position = if (relative) "fill" else "stack", width = 0.75) +
    ggplot2::scale_fill_manual(values = cfg$color_pals$nerve_subclusters) +
    {if (relative) ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1))} +
    ggplot2::labs(x = gsub("_", " ", outcome_col), y = if (relative) paste0(tools::toTitleCase(stat), " subcluster proportion") else paste0(tools::toTitleCase(stat), " STaN+ area (%)"), fill = "Subcluster") +
    cowplot::theme_cowplot()
  .save_subcluster_plot(p, filename, cfg, 5, 5); invisible(list(plot = p, summary = sum_df, data = long))
}
