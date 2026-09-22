#' nerve_clustering.R
#'
#' Two-stage k-means clustering pipeline used to i.) filter out imaging artifacts (non-neuronal cells expressing SYP/PGP9.5, RBCs, folded tissue, etc.), and ii.) identify and phenotype nerve objects
#' among all segmented multiplex-imaging objects, and to sub-type those
#' nerve objects into marker-based subclusters.
#'
#' @section Stage 1 (OG clustering):
#' k-means with k = 3 is run on all segmented objects passing
#' quality-control filtering, using all ten normalized marker intensity
#' channels. The resulting clusters are screened for a PGP9.5-enriched,
#' non-autofluorescence or DAPI-enriched signature; the cluster satisfying this criterion is
#' treated as the nerve-candidate population. Cluster identity is
#' determined by marker enrichment rather than a fixed numeric label,
#' since k-means cluster indices are arbitrary from run to run.
#'
#' @section Stage 2 (subclustering):
#' k-means with k = 7 is run on the nerve-candidate objects only, using
#' the same ten marker channels (scaled using the mean/SD of the full
#' object population from Stage 1, not re-scaled on the nerve-only
#' subset). One resulting subcluster represents a non-neuronal
#' false-positive population, identified by visual inspection of the
#' per-subcluster marker-enrichment and UMAP diagnostics (see
#' [plot_subcluster_enrichment()]/[plot_subcluster_umaps()]) rather than
#' by an automatic marker threshold, and is excluded from downstream nerve
#' quantification; the remaining six subclusters are retained.
#'
#' @section Input:
#' \describe{
#'   \item{object_count_df}{per-object marker intensity table (one row per
#'     segmented object), with a "Tissue Category", "Sample Name",
#'     "Category Region ID.obj_count_data", "Object ID.obj_count_data",
#'     and the ten marker "Mean (...)" columns. Images excluded during
#'     image-level quality review are assumed already removed from this
#'     table -- see [build_objs4pca()]'s "Note on input data" section.}
#'   \item{tissue_seg_df}{per-region tissue segmentation table.}
#' }
#'
#' @section Output:
#' A per-object data frame (`nerve_obj_table`) with `phenotype`
#' (`"subcluster_1".."subcluster_7"` when using the published labels, or
#' `"cluster_1".."cluster_<k>"` for a fresh from-scratch clustering run),
#' `cell_type` ("Nerve"/"Stroma"), `object_tag`, `region_tag`, and the
#' original marker/quantification columns.
#'
#' @section Reproducibility:
#' Exact reproduction of the published analyses uses archived per-object
#' Stage 2 cluster assignments from the original analysis, matched by
#' `object_tag`. Alternatively, Stage 2 k-means clustering can be performed
#' de novo. Because k-means solutions may vary across computational
#' environments, de novo clustering is intended for reanalysis rather than
#' exact reproduction of the published cluster assignments. Stage 1's
#' nerve-candidate cluster is identified dynamically by marker enrichment
#' rather than a fixed numeric label, so it is unaffected by this
#' distinction (see [run_og_kmeans()]).
NULL

library(dplyr)
library(tidyr)
library(rlang)
library(stringr)

## Parameters (marker channels, seed, k values, region filter, Stage 2
## settings, ...) are defined in config.R, which must be sourced first.


#' Compute region-level object density and identify excluded regions
#'
#' Aggregates object counts per region and joins them onto the full region
#' list from `tissue_seg_df`, so that regions with zero objects are still
#' represented. Used to identify regions excluded by the region-area/
#' density pre-filter (`region_filter_expr`) ahead of Stage 1 clustering.
#'
#' @param object_count_df Per-object marker intensity table (raw object-count data).
#' @param tissue_seg_df Per-region tissue segmentation table.
#' @param region_filter_expr An `rlang::expr()` evaluated against the
#'   region-level data frame (with `og_region_obj_density_per_mm2` and
#'   `` `Region Area (square microns).tissue_seg_data` `` in scope);
#'   regions satisfying this expression are *retained*, so the returned
#'   data frame is the *complement* (`!(!!region_filter_expr)`) -- i.e.
#'   the regions to omit.
#'
#' @return A data frame of regions that fail `region_filter_expr` (the
#'   regions to omit), with `region_tag`, object counts, and computed
#'   `og_region_obj_density_per_mm2`.
compute_regions_to_omit <- function(object_count_df, tissue_seg_df, region_filter_expr) {

  regions_plus_obj_density <- object_count_df %>%
    dplyr::select(`Sample Name`, `Tissue Category`,
                  `Category Region ID.obj_count_data`, `Object ID.obj_count_data`) %>%
    dplyr::filter(`Tissue Category` == "Stroma") %>%
    dplyr::left_join(
      tissue_seg_df %>%
        dplyr::select(`Sample Name`, `Tissue Category`,
                      `Region Area (square microns).tissue_seg_data`,
                      `Region ID.tissue_seg_data`),
      by = c("Sample Name", "Tissue Category",
             "Category Region ID.obj_count_data" = "Region ID.tissue_seg_data")
    ) %>%
    dplyr::mutate(region_tag = paste0(`Sample Name`, "_", `Tissue Category`, "_",
                                       `Category Region ID.obj_count_data`)) %>%
    dplyr::group_by(`Sample Name`, `Tissue Category`,
                    `Category Region ID.obj_count_data`, `region_tag`) %>%
    dplyr::summarize(
      og_n_total_obj_per_region = sum(!is.na(`Object ID.obj_count_data`)),
      `Region Area (square microns).tissue_seg_data` =
        dplyr::first(`Region Area (square microns).tissue_seg_data`),
      .groups = "drop") %>%
    dplyr::ungroup() %>%
    dplyr::right_join(
      (tissue_seg_df %>%
         dplyr::select(`Sample Name`, `Tissue Category`,
                       `Region Area (square microns).tissue_seg_data`,
                       `Region ID.tissue_seg_data`) %>%
         dplyr::mutate(region_tag = paste0(`Sample Name`, "_", `Tissue Category`, "_",
                                            `Region ID.tissue_seg_data`))),
      by = c("Sample Name", "Tissue Category",
             "Category Region ID.obj_count_data" = "Region ID.tissue_seg_data",
             "Region Area (square microns).tissue_seg_data", "region_tag")) %>%
    dplyr::mutate(
      og_n_total_obj_per_region = tidyr::replace_na(og_n_total_obj_per_region, 0),
      og_region_obj_density_per_mm2 =
        (og_n_total_obj_per_region / `Region Area (square microns).tissue_seg_data`) * 1e6)

  regions_plus_obj_density %>%
    dplyr::filter(!(!!region_filter_expr))
}


#' Stage 0: object-level filtering and tagging
#'
#' Tags each object/region with a stable identifier and restricts to the
#' target tissue category and non-omitted regions, ahead of Stage 1
#' clustering.
#'
#' @section Note on input data:
#' Images excluded during image-level quality review (omitted, redacted
#' due to imaging artifacts, or from non-target tissue types) are assumed
#' to already be removed from `object_count_df`/`tissue_seg_df` before
#' this function is called; no decision/category filtering is performed
#' here. The only pre-clustering filter applied is
#' [cfg$clustering$default_region_filter_expr] (region area > 2000 um^2). If you add new
#' raw data that hasn't been through the same image-level review, remove the
#' images excluded by that review before calling this function. The final
#' quantification cohort (`include_final_imaging_cohort`) is applied
#' downstream in `nerve_aggregation.R`, not at Stage 1.
#'
#' @param object_count_df Per-object marker intensity table (raw object-count data).
#' @param tissue_seg_df Per-region tissue segmentation table.
#' @param region_filter_expr An `rlang::expr()` of region-level filter
#'   conditions, passed to [compute_regions_to_omit()]. Defaults to
#'   [cfg$clustering$default_region_filter_expr].
#' @param quant_cols Character vector of the marker intensity columns to
#'   carry through to clustering. Defaults to [cfg$clustering$quant_cols].
#' @param tissue_keep Tissue category to retain (default `"Stroma"`).
#' @param id_col Column in `object_count_df` identifying the image
#'   (default `"Sample Name"`).
#'
#' @return A list with:
#' \describe{
#'   \item{objs4pca}{The filtered, tagged object table ready for Stage 1
#'     clustering.}
#'   \item{regions_to_omit}{Output of [compute_regions_to_omit()].}
#'   \item{quant_cols}{The `quant_cols` argument, passed through.}
#' }
build_objs4pca <- function(object_count_df,
                            tissue_seg_df,
                            region_filter_expr = cfg$clustering$default_region_filter_expr,
                            quant_cols = cfg$clustering$quant_cols,
                            tissue_keep = "Stroma",
                            id_col = "Sample Name") {

  qual_cols <- c("Sample Name", "Tissue Category", "object_tag", "TMA")
  keep_cols <- unique(c(quant_cols, qual_cols, "AF_to_PGP9.5+SYP_ratio.obj_count_data"))

  regions_to_omit <- compute_regions_to_omit(object_count_df, tissue_seg_df, region_filter_expr)

  # An `object_tag` column supplied in the input (e.g. from Table S4) is
  # used as-is; otherwise it is built from the source table's row names.
  has_tag <- "object_tag" %in% names(object_count_df)

  objs4pca <- object_count_df %>%
    dplyr::mutate(
      object_tag = if (has_tag) .data[["object_tag"]] else
        paste0(`Sample Name`, "_", `Tissue Category`, "_region",
               `Category Region ID.obj_count_data`, "_object",
               `Object ID.obj_count_data`, "_", row.names(.)),
      region_tag = paste0(`Sample Name`, "_", `Tissue Category`, "_",
                          `Category Region ID.obj_count_data`)
    ) %>%
    dplyr::filter(
      .data[["Tissue Category"]] == tissue_keep,
      !region_tag %in% regions_to_omit$region_tag
    ) %>%
    dplyr::mutate(TMA = paste0("TMA", stringr::str_split_i(.data[[id_col]], " ", 5))) %>%
    dplyr::select(dplyr::any_of(keep_cols)) %>%
    tidyr::drop_na(dplyr::any_of(quant_cols)) %>%
    as.data.frame()

  list(objs4pca = objs4pca, regions_to_omit = regions_to_omit, quant_cols = quant_cols)
}


#' Stage 1: k = 3 clustering of all filtered objects
#'
#' Runs k-means (k = `k`) on the scaled marker matrix of `objs4pca`, then
#' computes cluster-level marker enrichment (cluster mean minus population
#' mean) to identify which numeric cluster corresponds to the
#' nerve-candidate signature (PGP9.5-enriched, not autofluorescence-flagged).
#'
#' @param objs4pca Output of [build_objs4pca()]$objs4pca.
#' @param quant_cols Character vector of marker columns to cluster on.
#'   Defaults to [cfg$clustering$quant_cols].
#' @param k Number of clusters (default [cfg$clustering$k_og]).
#'
#' @return A list with:
#' \describe{
#'   \item{objs4pca_k}{`objs4pca` with cluster assignment columns added.}
#'   \item{scaled_objs4pca_k}{The scaled marker matrix (same rows/cluster
#'     labels), carried forward so Stage 2 subclustering reuses this
#'     population's scaling rather than re-scaling on the nerve-only subset.}
#'   \item{cluster_marker_scores}{Per-cluster, per-marker enrichment
#'     (`delta`), with `af_flag` and `possibly_nerves` annotations.}
#'   \item{kmeans}{The raw `stats::kmeans()` result.}
#'   \item{quant_cols, k}{Passed through.}
#' }
run_og_kmeans <- function(objs4pca, quant_cols = cfg$clustering$quant_cols, k = cfg$clustering$k_og) {

  scaled_data <- objs4pca
  scaled_data[, quant_cols] <- scale(objs4pca[, quant_cols])

  km <- stats::kmeans(scaled_data[, quant_cols], centers = k)

  objs4pca_k <- objs4pca %>%
    dplyr::mutate(!!paste0("preliminary_cluster_k", k) := km$cluster,
                  cluster = factor(km$cluster)) %>%
    as.data.frame()

  # The scaled marker matrix (with the same cluster labels attached) is
  # carried forward for subclustering, so that the nerve-only subset used
  # in Stage 2 is scaled using the full population's mean/SD rather than
  # being re-scaled on its own.
  scaled_objs4pca_k <- scaled_data %>%
    dplyr::mutate(!!paste0("preliminary_cluster_k", k) := km$cluster,
                  cluster = factor(km$cluster)) %>%
    as.data.frame()

  # Cluster-level marker enrichment (cluster mean minus population mean).
  cluster_means <- objs4pca_k %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(dplyr::across(dplyr::all_of(quant_cols), ~ mean(.x, na.rm = TRUE),
                                    .names = "mean_{.col}"), .groups = "drop")
  global_means <- objs4pca_k %>%
    dplyr::summarize(dplyr::across(dplyr::all_of(quant_cols), ~ mean(.x, na.rm = TRUE),
                                    .names = "global_{.col}"))

  cluster_marker_scores <- cluster_means %>%
    tidyr::pivot_longer(dplyr::starts_with("mean_"), names_to = "marker",
                        names_prefix = "mean_", values_to = "cluster_mean") %>%
    dplyr::left_join(
      global_means %>%
        tidyr::pivot_longer(dplyr::everything(), names_to = "marker",
                            names_prefix = "global_", values_to = "global_mean"),
      by = "marker") %>%
    dplyr::mutate(delta = cluster_mean - global_mean) %>%
    dplyr::arrange(cluster, dplyr::desc(delta)) %>%
    dplyr::mutate(
      marker_clean = marker |>
        gsub("Object | \\(.*?\\)| Mean.*|\\.obj_count_data|uto|luorescence", "", x = _) |>
        gsub("\\s+", " ", x = _)
    )

  # A cluster is flagged as autofluorescence-dominated if its autofluorescence
  # signal is substantially above the population mean; the nerve-candidate
  # cluster is defined as PGP9.5-enriched and not autofluorescence-flagged.
  cluster_marker_scores <- cluster_marker_scores %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(
      af_flag = any(marker_clean == "Af" & delta > 5),
      possibly_nerves = dplyr::case_when(
        marker_clean == "PGP9.5" & delta > 0 & !af_flag ~ "*",
        TRUE ~ NA_character_
      )) %>%
    dplyr::ungroup()

  list(objs4pca_k = objs4pca_k, scaled_objs4pca_k = scaled_objs4pca_k,
       cluster_marker_scores = as.data.frame(cluster_marker_scores),
       kmeans = km, quant_cols = quant_cols, k = k)
}

#' Identify the Stage 1 nerve-candidate cluster(s)
#'
#' @param cluster_marker_scores Output of [run_og_kmeans()]$cluster_marker_scores.
#'
#' @return Numeric vector of the Stage 1 cluster index/indices satisfying
#'   the PGP9.5-enrichment/non-autofluorescence nerve-candidate criteria.
select_nerve_og_cluster <- function(cluster_marker_scores) {
  cluster_marker_scores %>%
    dplyr::filter(af_flag == FALSE & possibly_nerves == "*") %>%
    dplyr::pull(cluster) %>%
    unique() %>%
    as.character() %>%
    as.numeric()
}

#' Annotate Stage 1 objects with a nerve-candidate flag and cluster label
#'
#' Adds two columns to `objs4pca_k`:
#' \describe{
#'   \item{candidate_nerve_obj}{`TRUE`/`FALSE` -- whether the object's OG
#'     cluster is one of `nerve_og_clusters` (output of
#'     [select_nerve_og_cluster()]).}
#'   \item{cluster_characteristic}{A human-readable `"<Marker>_high"` label
#'     for the object's OG cluster, based on whichever marker has the
#'     largest positive mean-vs-global enrichment `delta` for that cluster
#'     (e.g. `"PGP95_high"` for the nerve-candidate cluster, `"AF_high"`
#'     for an autofluorescence-dominated cluster). This describes the
#'     CLUSTER's dominant signature, not a per-object marker value.}
#'
#' @param objs4pca_k Output of [run_og_kmeans()]$objs4pca_k.
#' @param cluster_marker_scores Output of [run_og_kmeans()]$cluster_marker_scores.
#' @param nerve_og_clusters Output of [select_nerve_og_cluster()].
#' @param which_k_og Stage 1 `k` used (selects the
#'   `preliminary_cluster_k<which_k_og>` column). Defaults to [cfg$clustering$which_k_og].
#'
#' @return `objs4pca_k` with `candidate_nerve_obj` and
#'   `cluster_characteristic` columns added.
annotate_og_cluster_characteristics <- function(objs4pca_k, cluster_marker_scores,
                                                 nerve_og_clusters, which_k_og = cfg$clustering$which_k_og) {
  og_col <- paste0("preliminary_cluster_k", which_k_og)

  dominant <- cluster_marker_scores %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(delta, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      marker_label = gsub("\\.", "", marker_clean),
      marker_label = dplyr::if_else(marker_label == "Af", "AF", marker_label)
    ) %>%
    dplyr::transmute(cluster = as.character(cluster),
                      cluster_characteristic = paste0(marker_label, "_high"))

  objs4pca_k %>%
    dplyr::mutate(
      candidate_nerve_obj = .data[[og_col]] %in% nerve_og_clusters,
      cluster_characteristic = dominant$cluster_characteristic[
        match(as.character(.data[[og_col]]), dominant$cluster)]
    )
}


#' Stage 2: de novo phenotypic clustering of nerve-candidate objects
#'
#' Restricts the Stage 1 object population to nerve candidates and performs
#' k-means clustering across `k_min:k_max`. Marker values retain the scaling
#' calculated across the full Stage 1 object population rather than being
#' re-scaled within the nerve-candidate subset.
#'
#' @param objs4pca_k Output of [run_og_kmeans()]$objs4pca_k.
#' @param scaled_objs4pca_k Output of [run_og_kmeans()]$scaled_objs4pca_k.
#' @param nerve_og_clusters Numeric vector of Stage 1 cluster index/indices
#'   to retain (output of [select_nerve_og_cluster()]).
#' @param which_k_og The Stage 1 `k` used, so the corresponding
#'   `preliminary_cluster_k<which_k_og>` column can be read. Defaults
#'   to [cfg$clustering$which_k_og].
#' @param quant_cols Character vector of marker columns to cluster on.
#'   Defaults to [cfg$clustering$quant_cols].
#' @param k_min,k_max Range of subclustering `k` values to compute
#'   (default [cfg$clustering$k_sub_min]:[cfg$clustering$k_sub_max]).
#' @param nstart,iter.max `stats::kmeans()` parameters (default 50, 100).
#'
#' @return A list with:
#' \describe{
#'   \item{nerve_objs}{The nerve-candidate object table, with one
#'     `nerve_candidate_cluster_k<k>` column added per value of `k`.}
#'   \item{out_by_k}{Named list (`"k_<k>"`) of the raw `kmeans` result and
#'     column name for each `k`.}
#'   \item{which_k_og}{Passed through.}
#'   \item{scaled_nerve}{The scaled marker matrix restricted to
#'     nerve-candidate objects.}
#' }
subcluster_nerve_objects <- function(objs4pca_k,
                                     scaled_objs4pca_k,
                                     nerve_og_clusters,
                                     which_k_og = cfg$clustering$which_k_og,
                                     quant_cols = cfg$clustering$quant_cols,
                                     k_min = cfg$clustering$k_sub_min,
                                     k_max = cfg$clustering$k_sub_max,
                                     nstart = cfg$clustering$nstart_denovo,
                                     iter.max = cfg$clustering$iter_max_denovo) {

  og_cluster_col <- paste0("preliminary_cluster_k", which_k_og)

  nerve_objs <- objs4pca_k %>%
    dplyr::filter(.data[[og_cluster_col]] %in% nerve_og_clusters)

  # Retain marker scaling calculated across the full Stage 1 population.
  scaled_nerve <- scaled_objs4pca_k %>%
    dplyr::filter(object_tag %in% nerve_objs$object_tag)

  out_by_k <- list()
  for (k in k_min:k_max) {
    km <- stats::kmeans(scaled_nerve[, quant_cols], centers = k,
                         nstart = nstart, iter.max = iter.max)

    col_name <- paste0("nerve_candidate_cluster_k", k)
    nerve_objs[[col_name]] <- factor(km$cluster)

    out_by_k[[paste0("k_", k)]] <- list(kmeans = km, colname = col_name)
  }

  list(nerve_objs = nerve_objs, out_by_k = out_by_k, which_k_og = which_k_og,
       scaled_nerve = scaled_nerve)
}

#' Validate archived publication assignments before use
#'
#' Checks that `object_tag` is present in both tables, unique in
#' `publication_reference`, and covers every nerve-candidate object.
#'
#' @param nerve_candidates Data frame of nerve-candidate objects, with an
#'   `object_tag` column.
#' @param publication_reference Data frame of archived per-object cluster
#'   assignments, with an `object_tag` column.
#' @param context Character string identifying the calling context, used
#'   in error messages.
#' @keywords internal
.validate_publication_reference <- function(nerve_candidates, publication_reference, context) {
  if (!"object_tag" %in% names(nerve_candidates)) {
    stop(context, ": nerve-candidate object table has no `object_tag` column.")
  }
  if (!"object_tag" %in% names(publication_reference)) {
    stop(context, ": publication_reference has no `object_tag` column.")
  }
  if (anyDuplicated(publication_reference$object_tag)) {
    stop(context, ": publication_reference$object_tag must be unique.")
  }
  missing <- setdiff(nerve_candidates$object_tag, publication_reference$object_tag)
  if (length(missing) > 0) {
    stop(context, ": ", length(missing), " of ", nrow(nerve_candidates),
         " nerve-candidate objects have no archived assignment in publication_reference ",
         "(first missing: ", missing[1], ").")
  }
}

#' Full pipeline: object filtering, Stage 1/2 clustering, and phenotyping
#'
#' Chains [build_objs4pca()], [run_og_kmeans()],
#' [select_nerve_og_cluster()], and [subcluster_nerve_objects()], then
#' derives the final per-object `phenotype`/`cell_type` columns.
#'
#' @param object_count_df Per-object marker intensity table (raw object-count
#'   data), already restricted to the reviewed image set -- see
#'   [build_objs4pca()]'s "Note on input data" section for why no
#'   decision/category filtering happens in this pipeline.
#' @param tissue_seg_df Per-region tissue segmentation table.
#' @param region_filter_expr An `rlang::expr()` of region-level filter
#'   conditions. Defaults to [cfg$clustering$default_region_filter_expr].
#' @param quant_cols Character vector of marker columns to cluster on.
#'   Defaults to [cfg$clustering$quant_cols].
#' @param seed Random seed for `set.seed()` before Stage 1 clustering.
#'   Defaults to [cfg$clustering$seed].
#' @param k_og Stage 1 number of clusters. Defaults to [cfg$clustering$k_og].
#' @param which_k_og Stage 1 `k` used to select nerve-candidate objects
#'   for Stage 2. Defaults to [cfg$clustering$which_k_og].
#' @param k_sub_min,k_sub_max Range of Stage 2 subclustering `k` values to
#'   compute. Default [cfg$clustering$k_sub_min]/[cfg$clustering$k_sub_max].
#' @param k_sub_final Which Stage 2 `k` to use for the final
#'   `phenotype`/`cell_type` columns. Defaults to [cfg$clustering$k_sub_max]. Must lie
#'   within `k_sub_min:k_sub_max`.
#' @param use_publication_cluster_labels Logical. If `TRUE`, use archived
#'   per-object Stage 2 cluster assignments from the original analysis for
#'   exact reproduction of published results. If `FALSE`, perform Stage 2
#'   k-means clustering de novo. Default `TRUE`.
#' @param omit_nonneuronal_cluster_3_from_publication Logical. Applies only
#'   when `use_publication_cluster_labels = FALSE` (when `TRUE`, `cell_type`
#'   already comes directly from `publication_reference`). If `TRUE`
#'   (default), `cell_type` is set to `"Stroma"` for the objects identified
#'   as the non-neuronal subcluster in the published analysis, matched by
#'   `object_tag` from `publication_reference`, independent of the de novo
#'   cluster labels. If `FALSE`, `cell_type` is instead determined by
#'   `nerve_subclusters_final` applied to the de novo labels.
#' @param publication_reference Data frame containing archived per-object
#'   cluster assignments from the original analysis. Assignments are
#'   matched to nerve objects by `object_tag`. Required when
#'   `use_publication_cluster_labels` or
#'   `omit_nonneuronal_cluster_3_from_publication` is `TRUE`.
#' @param nerve_subclusters_final Numeric vector of subcluster indices (a
#'   subset of `1:k_sub_final`) retained as `cell_type == "Nerve"`. Used
#'   only when both `use_publication_cluster_labels` and
#'   `omit_nonneuronal_cluster_3_from_publication` are `FALSE`, in which
#'   case it must be supplied explicitly based on this run's own Stage 2
#'   diagnostics (see [plot_subcluster_enrichment()]/
#'   [plot_subcluster_umaps()]). [cfg$clustering$nerve_subclusters_final] gives the
#'   published k = 7 assignment and applies only to that clustering, not to
#'   a de novo run.
#' @param nstart,iter.max Passed to [subcluster_nerve_objects()] when
#'   `use_publication_cluster_labels = FALSE`. Default 50, 100.
#'
#' @return A list with:
#' \describe{
#'   \item{objs4pca_og_clustered}{Output of [run_og_kmeans()]$objs4pca_k,
#'     with `candidate_nerve_obj` and `cluster_characteristic` added via
#'     [annotate_og_cluster_characteristics()].}
#'   \item{scaled_objs4pca_og_clustered}{Output of [run_og_kmeans()]$scaled_objs4pca_k.}
#'   \item{og_cluster_marker_scores}{Output of [run_og_kmeans()]$cluster_marker_scores.}
#'   \item{nerve_og_clusters}{Output of [select_nerve_og_cluster()].}
#'   \item{subcluster_result}{Output of [subcluster_nerve_objects()] when
#'     `use_publication_cluster_labels = FALSE`; `NULL` when `TRUE`, since
#'     de novo Stage 2 clustering is not performed in that path.}
#'   \item{nerve_obj_table}{Final per-object table with `subcluster_num`,
#'     `phenotype`, and `cell_type` columns -- the primary input to
#'     `nerve_aggregation.R`.}
#'   \item{quant_cols}{Passed through.}
#'   \item{regions_to_omit}{Output of [build_objs4pca()]$regions_to_omit;
#'     pass this through to `nerve_aggregation.R` so the same
#'     regions are excluded at the final aggregation step.}
#' }
#'
#' @examples
#' \dontrun{
#' # Reproduce the published phenotyping exactly (default):
#' result <- run_nerve_object_kmeans_pipeline(
#'   object_count_df       = obj_count_data,
#'   tissue_seg_df         = tissue_seg_data,
#'   publication_reference = readRDS(cfg$paths$publication_stage2_reference)
#' )
#' nerve_obj_table <- result$nerve_obj_table
#'
#' # De novo Stage 2 clustering (e.g. exploring a different k), independent
#' # of the published assignments:
#' result_denovo <- run_nerve_object_kmeans_pipeline(
#'   object_count_df = obj_count_data,
#'   tissue_seg_df = tissue_seg_data,
#'   k_sub_min = 3, k_sub_max = 5, k_sub_final = 5,
#'   use_publication_cluster_labels = FALSE,
#'   omit_nonneuronal_cluster_3_from_publication = FALSE,
#'   nerve_subclusters_final = c(1, 2, 4, 5)
#' )
#' }
run_nerve_object_kmeans_pipeline <- function(object_count_df,
                                             tissue_seg_df,
                                             region_filter_expr = cfg$clustering$default_region_filter_expr,
                                             quant_cols = cfg$clustering$quant_cols,
                                             seed = cfg$clustering$seed,
                                             k_og = cfg$clustering$k_og,
                                             which_k_og = cfg$clustering$which_k_og,
                                             k_sub_min = cfg$clustering$k_sub_min,
                                             k_sub_max = cfg$clustering$k_sub_max,
                                             k_sub_final = cfg$clustering$k_sub_max,
                                             use_publication_cluster_labels = cfg$clustering$use_publication_cluster_labels,
                                             omit_nonneuronal_cluster_3_from_publication = cfg$clustering$omit_nonneuronal_cluster_3_from_publication,
                                             publication_reference = NULL,
                                             nerve_subclusters_final = NULL,
                                             nstart = cfg$clustering$nstart_denovo,
                                             iter.max = cfg$clustering$iter_max_denovo) {

  stopifnot(
    "k_sub_final must be between k_sub_min and k_sub_max" =
      k_sub_final >= k_sub_min && k_sub_final <= k_sub_max
  )

  needs_reference <- use_publication_cluster_labels || omit_nonneuronal_cluster_3_from_publication
  if (needs_reference && is.null(publication_reference)) {
    stop("publication_reference must be supplied (a data frame with `object_tag`, ",
         "`subcluster_num`, `phenotype`, `cell_type` columns) when ",
         "use_publication_cluster_labels or omit_nonneuronal_cluster_3_from_publication ",
         "is TRUE -- load it with readRDS(cfg$paths$publication_stage2_reference) and pass it ",
         "in, or set both to FALSE for de novo clustering.")
  }
  if (!use_publication_cluster_labels && !omit_nonneuronal_cluster_3_from_publication &&
      is.null(nerve_subclusters_final)) {
    stop("nerve_subclusters_final must be supplied when both ",
         "use_publication_cluster_labels and omit_nonneuronal_cluster_3_from_publication ",
         "are FALSE -- inspect this run's Stage 2 diagnostics (see ",
         "plot_subcluster_enrichment()/plot_subcluster_umaps()) and supply the indices to ",
         "retain as cell_type == \"Nerve\".")
  }
  if (!is.null(nerve_subclusters_final) && !all(nerve_subclusters_final %in% seq_len(k_sub_final))) {
    stop("nerve_subclusters_final must be a subset of 1:k_sub_final.")
  }

  set.seed(seed)

  step1 <- build_objs4pca(object_count_df, tissue_seg_df,
                          region_filter_expr = region_filter_expr,
                          quant_cols = quant_cols)

  og <- run_og_kmeans(step1$objs4pca, quant_cols = quant_cols, k = k_og)

  nerve_og_clusters <- select_nerve_og_cluster(og$cluster_marker_scores)
  if (length(nerve_og_clusters) == 0) {
    stop("No Stage 1 cluster satisfied the PGP9.5-enrichment/autofluorescence ",
         "nerve-candidate criteria -- check region_filter_expr and input data.")
  }
  message("Nerve-candidate cluster(s): ", paste(nerve_og_clusters, collapse = ", "))

  og_objs4pca_k_annotated <- annotate_og_cluster_characteristics(
    og$objs4pca_k, og$cluster_marker_scores, nerve_og_clusters, which_k_og = which_k_og)

  if (use_publication_cluster_labels) {
    # Publication reproduction: apply archived Stage 2 assignments by
    # object_tag. De novo Stage 2 clustering is not performed in this path.
    og_cluster_col <- paste0("preliminary_cluster_k", which_k_og)
    nerve_obj_table <- og$objs4pca_k %>%
      dplyr::filter(.data[[og_cluster_col]] %in% nerve_og_clusters)

    .validate_publication_reference(nerve_obj_table, publication_reference,
                                     context = "use_publication_cluster_labels = TRUE")

    # Match archived publication assignments by stable object identifier.
    matched <- publication_reference[
      match(nerve_obj_table$object_tag, publication_reference$object_tag),
      c("subcluster_num", "phenotype", "cell_type")]
    nerve_obj_table$subcluster_num <- matched$subcluster_num
    nerve_obj_table$phenotype <- matched$phenotype
    nerve_obj_table$cell_type <- matched$cell_type

    sub <- NULL
  } else {
    # De novo clustering: perform Stage 2 k-means and derive phenotype/
    # cell_type from this run's own labels.
    sub <- subcluster_nerve_objects(og$objs4pca_k, og$scaled_objs4pca_k, nerve_og_clusters,
                                    which_k_og = which_k_og, quant_cols = quant_cols,
                                    k_min = k_sub_min, k_max = k_sub_max,
                                    nstart = nstart, iter.max = iter.max)

    final_col <- paste0("nerve_candidate_cluster_k", k_sub_final)
    nerve_obj_table <- sub$nerve_objs
    nerve_obj_table$subcluster_num <- as.integer(as.character(nerve_obj_table[[final_col]]))
    nerve_obj_table$phenotype <- paste0("cluster_", nerve_obj_table$subcluster_num)

    if (omit_nonneuronal_cluster_3_from_publication) {
      .validate_publication_reference(nerve_obj_table, publication_reference,
                                       context = "omit_nonneuronal_cluster_3_from_publication = TRUE")
      nerve_obj_table$cell_type <- publication_reference$cell_type[
        match(nerve_obj_table$object_tag, publication_reference$object_tag)]
    } else {
      nerve_obj_table$cell_type <- dplyr::if_else(
        nerve_obj_table$subcluster_num %in% nerve_subclusters_final, "Nerve", "Stroma")
    }
  }

  list(
    objs4pca_og_clustered = og_objs4pca_k_annotated,
    scaled_objs4pca_og_clustered = og$scaled_objs4pca_k,
    og_cluster_marker_scores = og$cluster_marker_scores,
    nerve_og_clusters = nerve_og_clusters,
    subcluster_result = sub,
    nerve_obj_table = nerve_obj_table,
    quant_cols = quant_cols,
    # Regions excluded by the region-area/density pre-filter; pass this
    # through to nerve_aggregation.R so the same regions are excluded
    # at the final aggregation step.
    regions_to_omit = step1$regions_to_omit
  )
}

#' Diagnostic plots
#'
#' For each clustering stage: PCA (`factoextra::fviz_cluster`) and UMAP
#' embeddings colored by cluster assignment, a marker-enrichment bar
#' chart, and a UMAP grid colored by each individual marker. The UMAP
#' embedding is computed only for visualization and has no bearing on
#' cluster assignment, which is determined by k-means on the (scaled)
#' marker matrix directly.
NULL

#' Marker display order for the Stage 1 enrichment bar chart.
MARKERS_ORDER <- c("PGP9.5", "SYP", "SP", "TH", "VACHT", "VGLUT1", "CD31", "DAPI", "panCK", "Af")

#' Nerve-subcluster colors, matching `cfg$color_pals$nerve_subclusters` in
#' `config.R` (used for the Stage 2 PCA/UMAP/enrichment plots
#' so subcluster colors are consistent with the rest of the paper's
#' figures).
NERVE_SUBCLUSTER_COLORS <- c(
  "1" = "#C4B8F2", "2" = "#F7C9D5", "3" = "#5D46F2", "4" = "#96E082",
  "5" = "#F2178E", "6" = "#FEC83E", "7" = "#F28700"
)

#' Plot labels for the k = 7 Stage 2 subclusters, keyed by subcluster number
#' (matches `cluster_labels` in `config.R`). Used for the
#' `cluster_phenotype` column of the supplementary object tables.
NERVE_SUBCLUSTER_LABELS <- c(
  "1" = "Cluster 1 (VGLUT1+/SP+)",
  "2" = "Cluster 2 (TH-mid)",
  "3" = "Cluster 3 (SYP+ non-neuronal)",
  "4" = "Cluster 4 (VAChT+)",
  "5" = "Cluster 5 (TH-high)",
  "6" = "Cluster 6 (SYP-high)",
  "7" = "Cluster 7 (Indeterminate/Other)"
)

#' Strip a marker column name down to its short display name
#'
#' @param x Character vector of full marker column names (e.g.
#'   `` "Object PGP9.5 (Opal 480) Mean (...).obj_count_data" ``).
#'
#' @return Character vector of short marker names (e.g. `"PGP9.5"`).
#' @keywords internal
.marker_shortname <- function(x) {
  gsub(" \\(.*", "", gsub("Object |Mean |\\.obj_count_data|uto|luorescence| \\(.*", "", x))
}

#' Compute a 2-D UMAP embedding, for visualization only
#'
#' @param scaled_df Data frame containing the scaled marker columns.
#' @param quant_cols Character vector of marker columns to embed.
#' @param seed Random seed for `set.seed()` before `uwot::umap()`.
#'   Defaults to [cfg$clustering$seed].
#' @param n_neighbors,min_dist `uwot::umap()` parameters.
#'
#' @return `scaled_df` with `UMAP1`/`UMAP2` columns added.
compute_umap_embedding <- function(scaled_df, quant_cols, seed = cfg$clustering$seed,
                                    n_neighbors = 15, min_dist = 0.0125) {
  set.seed(seed)
  um <- uwot::umap(scaled_df[, quant_cols], n_neighbors = n_neighbors,
                    min_dist = min_dist, metric = "euclidean")
  as.data.frame(um) |>
    dplyr::rename(UMAP1 = V1, UMAP2 = V2) |>
    dplyr::bind_cols(scaled_df)
}

#' Stage 1 diagnostic plots: PCA, UMAP, enrichment, and per-marker UMAPs
#'
#' @param objs4pca_k Output of [run_og_kmeans()]$objs4pca_k.
#' @param scaled_objs4pca_k Output of [run_og_kmeans()]$scaled_objs4pca_k.
#' @param cluster_marker_scores Output of [run_og_kmeans()]$cluster_marker_scores.
#' @param quant_cols Character vector of marker columns. Defaults to [cfg$clustering$quant_cols].
#' @param n_objs_to_plot Number of objects to subsample for plotting
#'   (default 6666).
#' @param markers_order Marker display order for the enrichment bar chart.
#'   Defaults to [MARKERS_ORDER].
#' @param seed Random seed for subsampling and UMAP. Defaults to [cfg$clustering$seed].
#'
#' @return A list with:
#' \describe{
#'   \item{pca_umap_enrichment}{A combined PCA + UMAP + enrichment-bar-chart
#'     `patchwork` plot.}
#'   \item{marker_umaps}{A `patchwork` grid of per-marker UMAPs.}
#' }
plot_og_clustering_diagnostics <- function(objs4pca_k, scaled_objs4pca_k,
                                            cluster_marker_scores,
                                            quant_cols = cfg$clustering$quant_cols,
                                            n_objs_to_plot = 6666,
                                            markers_order = MARKERS_ORDER,
                                            seed = cfg$clustering$seed) {
  set.seed(seed)
  objs2plot <- dplyr::slice_sample(scaled_objs4pca_k, n = min(n_objs_to_plot, nrow(scaled_objs4pca_k)))
  umap_df <- compute_umap_embedding(objs2plot, quant_cols, seed = seed) |>
    dplyr::mutate(cluster = factor(cluster))

  pca_plot <- factoextra::fviz_cluster(
    list(data = objs2plot[, quant_cols], cluster = objs2plot$cluster),
    labelsize = 0, alpha = 0.2, pointsize = 1.4
  ) +
    ggplot2::labs(title = "PCA of marker expression") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(aspect.ratio = 1, panel.grid = ggplot2::element_blank()) +
    ggplot2::coord_equal()

  umap_plot <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color = cluster, shape = cluster)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.8, na.rm = TRUE) +
    ggplot2::labs(title = "UMAP of marker expression") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(aspect.ratio = 1, panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

  enrichment_plot <- ggplot2::ggplot(
    cluster_marker_scores,
    ggplot2::aes(x = factor(marker_clean, levels = markers_order), y = delta,
                 fill = cluster, label = possibly_nerves)
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0, color = "grey20") +
    ggplot2::facet_wrap(~cluster, scales = "free_y", ncol = 1) +
    ggplot2::labs(title = "Enriched markers per cluster",
                  caption = "* = potential nerve cluster",
                  y = "Expression delta (cluster - global)") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   legend.position = "none")

  marker_umaps <- lapply(quant_cols, function(m) {
    q <- stats::quantile(umap_df[[m]], c(0.05, 0.95), na.rm = TRUE)
    ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color = .data[[m]])) +
      ggplot2::geom_point(alpha = 0.3, size = 0.25, na.rm = TRUE) +
      ggplot2::scale_color_gradientn(name = .marker_shortname(m), limits = q,
                                      colours = cetcolor::cet_pal(9, name = "l8"),
                                      oob = scales::squish) +
      ggplot2::theme_bw(base_size = 7) +
      ggplot2::theme(aspect.ratio = 1, panel.grid = ggplot2::element_blank(),
                     axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                     plot.title = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
  })
  names(marker_umaps) <- vapply(quant_cols, .marker_shortname, character(1))

  list(
    pca_umap_enrichment = (pca_plot / umap_plot) | enrichment_plot,
    marker_umaps = patchwork::wrap_plots(marker_umaps, ncol = 3) +
      patchwork::plot_annotation(title = "UMAP of marker expression (scaled)")
  )
}

#' Markers shown on the subcluster enrichment/UMAP figures (a subset of
#' the 10 clustering markers -- matches
#' `cfg$nerve_subclusters$markers_keep`/`umap_plot_markers` in
#' `config.R`; CD31/DAPI/panCK/Af are used for clustering but
#' not shown on these particular figures).
SUBCLUSTER_ENRICHMENT_MARKERS <- c("PGP9.5", "SYP", "SP", "TH", "VACHT", "VGLUT1")
SUBCLUSTER_UMAP_MARKERS <- c("PGP9.5", "SYP", "TH", "VACHT", "SP", "VGLUT1")

#' Cap and rescale a marker to \[0, 1\] for UMAP color scales
#'
#' Caps `x` at its `q`-th quantile, then min-max rescales to \[0, 1\].
#' Matches `cfg$nerve_subclusters$marker_cap_quantile` usage
#' (`.cap_minmax()` in `R/nerve_cluster_analysis.R`) -- this is the
#' "Scaled MPI" shown on the marker-expression UMAPs, distinct from the
#' z-scored matrix used for clustering itself.
#'
#' @param x Numeric vector.
#' @param q Quantile (0-1) at which to cap `x` before rescaling.
#'
#' @return Numeric vector rescaled to \[0, 1\] (all zero if `x` has no
#'   finite range).
cap_minmax <- function(x, q) {
  x <- pmin(x, stats::quantile(x, q, na.rm = TRUE))
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

#' Stage 2 marker-enrichment bar chart
#'
#' @param nerve_objs Nerve-candidate object table with a subcluster
#'   assignment column (e.g. [subcluster_nerve_objects()]$nerve_objs).
#' @param cluster_col Name of the subcluster assignment column in
#'   `nerve_objs`.
#' @param quant_cols Character vector of marker columns. Defaults to [cfg$clustering$quant_cols].
#' @param markers_keep Markers to display. Defaults to
#'   [SUBCLUSTER_ENRICHMENT_MARKERS].
#' @param subcluster_colors Named color vector keyed by subcluster.
#'   Defaults to [NERVE_SUBCLUSTER_COLORS].
#'
#' @return A `ggplot2` object.
plot_subcluster_enrichment <- function(nerve_objs, cluster_col,
                                        quant_cols = cfg$clustering$quant_cols,
                                        markers_keep = SUBCLUSTER_ENRICHMENT_MARKERS,
                                        subcluster_colors = NERVE_SUBCLUSTER_COLORS) {
  cluster_means <- nerve_objs |>
    dplyr::group_by(.data[[cluster_col]]) |>
    dplyr::summarize(dplyr::across(dplyr::all_of(quant_cols), ~ mean(.x, na.rm = TRUE),
                                    .names = "mean_{.col}"), .groups = "drop")
  global_means <- nerve_objs |>
    dplyr::summarize(dplyr::across(dplyr::all_of(quant_cols), ~ mean(.x, na.rm = TRUE),
                                    .names = "global_{.col}"))
  cluster_marker_scores <- cluster_means |>
    tidyr::pivot_longer(dplyr::starts_with("mean_"), names_to = "marker",
                        names_prefix = "mean_", values_to = "cluster_mean") |>
    dplyr::left_join(
      global_means |>
        tidyr::pivot_longer(dplyr::everything(), names_to = "marker",
                            names_prefix = "global_", values_to = "global_mean"),
      by = "marker") |>
    dplyr::mutate(delta = cluster_mean - global_mean, marker_clean = .marker_shortname(marker)) |>
    dplyr::filter(marker_clean %in% markers_keep) |>
    dplyr::mutate(marker_clean = factor(marker_clean, levels = markers_keep)) |>
    dplyr::rename(cluster = !!cluster_col)

  ggplot2::ggplot(cluster_marker_scores, ggplot2::aes(x = marker_clean, y = delta, fill = cluster)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0, color = "grey20") +
    ggplot2::facet_wrap(~cluster, scales = "free_y", ncol = 1) +
    ggplot2::scale_fill_manual(values = subcluster_colors) +
    ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    ggplot2::labs(title = "Enriched markers per cluster", x = NULL,
                  y = "Expression delta (cluster - global)") +
    ggplot2::coord_cartesian(clip = "off") +
    cowplot::theme_cowplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1.1),
                   panel.spacing.y = grid::unit(0.08, "cm")) +
    ggplot2::guides(fill = "none")
}

#' Stage 2 UMAP: cluster assignment plus per-marker expression grid
#'
#' @param nerve_objs Nerve-candidate object table with a subcluster
#'   assignment column.
#' @param scaled_nerve Scaled marker matrix for the same objects (e.g.
#'   [subcluster_nerve_objects()]$scaled_nerve).
#' @param cluster_col Name of the subcluster assignment column in
#'   `nerve_objs`.
#' @param quant_cols Character vector of marker columns. Defaults to [cfg$clustering$quant_cols].
#' @param umap_plot_markers Markers to display on the per-marker UMAP
#'   grid. Defaults to [SUBCLUSTER_UMAP_MARKERS].
#' @param marker_cap_quantile Quantile passed to [cap_minmax()] for the
#'   "Scaled MPI" color scale (default 0.98).
#' @param subcluster_colors Named color vector keyed by subcluster.
#'   Defaults to [NERVE_SUBCLUSTER_COLORS].
#' @param seed Random seed for UMAP. Defaults to [cfg$clustering$seed].
#'
#' @return A list with:
#' \describe{
#'   \item{cluster_plot}{UMAP colored by subcluster assignment.}
#'   \item{marker_plots}{List of per-marker UMAP `ggplot2` objects.}
#'   \item{combined_plot}{`cluster_plot` and `marker_plots` combined via `patchwork`.}
#' }
plot_subcluster_umaps <- function(nerve_objs, scaled_nerve, cluster_col,
                                   quant_cols = cfg$clustering$quant_cols,
                                   umap_plot_markers = SUBCLUSTER_UMAP_MARKERS,
                                   marker_cap_quantile = 0.98,
                                   subcluster_colors = NERVE_SUBCLUSTER_COLORS,
                                   seed = cfg$clustering$seed) {
  scaled_nerve$cluster <- factor(nerve_objs[[cluster_col]][match(scaled_nerve$object_tag, nerve_objs$object_tag)])
  umap_df <- compute_umap_embedding(scaled_nerve, quant_cols, seed = seed)

  marker_cols <- quant_cols[vapply(quant_cols, .marker_shortname, character(1)) %in% umap_plot_markers]
  names(marker_cols) <- vapply(marker_cols, .marker_shortname, character(1))
  for (nm in names(marker_cols)) {
    umap_df[[nm]] <- cap_minmax(umap_df[[marker_cols[nm]]], marker_cap_quantile)
  }

  cluster_plot <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color = cluster)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.9) +
    ggplot2::scale_color_manual(values = subcluster_colors) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "UMAP of nerve objects", color = "Subcluster") +
    cowplot::theme_cowplot() +
    ggplot2::theme(aspect.ratio = 1, axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))

  marker_plots <- lapply(names(marker_cols), function(nm) {
    ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color = .data[[nm]])) +
      ggplot2::geom_point(size = 0.2, alpha = 0.8) +
      ggplot2::scale_color_gradientn(name = "Scaled\nMPI", limits = c(0, 1),
                                      values = scales::rescale(c(0.15, 1)),
                                      colours = cetcolor::cet_pal(9, name = "l8"),
                                      oob = scales::squish) +
      ggplot2::coord_equal() +
      ggplot2::labs(title = nm) +
      cowplot::theme_cowplot() +
      ggplot2::theme(aspect.ratio = 1, axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
  })

  list(
    cluster_plot = cluster_plot,
    marker_plots = marker_plots,
    combined_plot = (patchwork::wrap_plots(cluster_plot,
                       patchwork::wrap_plots(marker_plots, ncol = 2, guides = "collect"),
                       ncol = 2) & ggplot2::theme(legend.position = "right"))
  )
}

#' Stage 2 UMAP for the subcluster figure (Fig. 6c-g)
#'
#' Generates the subcluster UMAP shown in Fig. 6c-g. It differs from
#' [plot_subcluster_umaps()] (a general diagnostic that embeds a single
#' population) in the following ways:
#'
#' \enumerate{
#'   \item The UMAP \strong{embedding} is computed from the Stage 1 scaled
#'     object table (`scaled_objs4pca_k`, restricted to nerve objects) for
#'     the images in `og_im3s`, using the raw z-scored marker values (no
#'     percentile cap or rescaling).
#'   \item The \strong{point color} (marker expression) is taken from a
#'     separate table built from `nerve_objs` for the images in
#'     `final_im3s` (the final patient-level cohort), with each marker
#'     independently capped at its `marker_cap_quantile` and min-max
#'     rescaled to \[0, 1\] via [cap_minmax()]. `final_im3s` and `og_im3s`
#'     are different image lists and both are required.
#'   \item The figure uses `seed = 444` and `repulsion_strength = 2`,
#'     which [compute_umap_embedding()] does not expose.
#' }
#'
#' Only points in `final_im3s` are plotted; the wider `og_im3s` population
#' contributes to the embedding geometry but is not itself shown.
#'
#' @param scaled_objs4pca_k Full OG-scaled Stage-1 candidate table, e.g.
#'   [run_og_kmeans()]$scaled_objs4pca_k -- NOT restricted to nerve objects.
#' @param nerve_objs Nerve-candidate object table with the FINAL subcluster
#'   assignment column, e.g. [subcluster_nerve_objects()]$nerve_objs.
#' @param cluster_col Name of the subcluster assignment column in `nerve_objs`.
#' @param final_im3s Character vector of `Sample Name` values from the FINAL
#'   patient-level nerve-quant table's `consolidated_im3s` column.
#' @param og_im3s Character vector of `Sample Name` values for the images
#'   used to compute the UMAP embedding (a different list from
#'   `final_im3s`).
#' @param quant_cols All ten marker columns. Defaults to [cfg$clustering$quant_cols].
#' @param umap_plot_markers The six markers shown for expression coloring.
#'   Defaults to [SUBCLUSTER_UMAP_MARKERS].
#' @param seed,n_neighbors,min_dist,repulsion_strength `uwot::umap()`
#'   parameters. Defaults (444, 15, 0.0125, 2) match the published figure.
#' @param marker_cap_quantile Quantile passed to [cap_minmax()] for the
#'   expression-coloring rescale (default 0.98).
#' @param subcluster_colors Named color vector keyed by subcluster.
#'   Defaults to [NERVE_SUBCLUSTER_COLORS].
#'
#' @return A list with `embedding_population` (`d` -- the `og_im3s`-filtered,
#'   ALL-subcluster population, with raw OG-scaled marker values, that was
#'   actually fed into `uwot::umap()` to compute the embedding), `umap_df`
#'   (that same population with `UMAP1`/`UMAP2` attached, all `og_im3s`
#'   objects), `umap_df_filt` (the `final_im3s` subset actually plotted,
#'   with recolored marker values), `cluster_plot`, `marker_plots`, and
#'   `combined_plot`.
plot_subcluster_umap_figure <- function(scaled_objs4pca_k, nerve_objs, cluster_col,
                                         final_im3s, og_im3s,
                                         quant_cols = cfg$clustering$quant_cols,
                                         umap_plot_markers = SUBCLUSTER_UMAP_MARKERS,
                                         seed = 444, n_neighbors = 15, min_dist = 0.0125,
                                         repulsion_strength = 2, marker_cap_quantile = 0.98,
                                         subcluster_colors = NERVE_SUBCLUSTER_COLORS) {
  marker_cols <- quant_cols
  names(marker_cols) <- vapply(quant_cols, .marker_shortname, character(1))
  nerve_tags <- nerve_objs$object_tag

  # embedding population: real nerve objects, og_im3s images, raw OG-scaled values
  d <- scaled_objs4pca_k[scaled_objs4pca_k$object_tag %in% nerve_tags &
                           scaled_objs4pca_k[["Sample Name"]] %in% og_im3s, ]
  for (nm in names(marker_cols)) d[[nm]] <- d[[marker_cols[[nm]]]]
  d[[cluster_col]] <- nerve_objs[[cluster_col]][match(d$object_tag, nerve_objs$object_tag)]

  # recolor population: final_im3s images, independently capped + rescaled
  d_ <- nerve_objs[nerve_objs[["Sample Name"]] %in% final_im3s, ]
  for (nm in names(marker_cols)) d_[[nm]] <- cap_minmax(d_[[marker_cols[[nm]]]], marker_cap_quantile)

  set.seed(seed)
  umap_out <- uwot::umap(d[, names(marker_cols)], n_neighbors = n_neighbors, min_dist = min_dist,
                          repulsion_strength = repulsion_strength, metric = "euclidean")

  umap_df <- as.data.frame(umap_out)
  names(umap_df) <- c("UMAP1", "UMAP2")
  umap_df[[cluster_col]] <- d[[cluster_col]]
  umap_df$object_tag <- d$object_tag
  umap_df[["Sample Name"]] <- d[["Sample Name"]]
  for (nm in umap_plot_markers) umap_df[[nm]] <- d[[nm]]

  umap_df_filt <- umap_df[umap_df[["Sample Name"]] %in% final_im3s, ]
  umap_df_filt[umap_plot_markers] <- NULL
  umap_df_filt <- merge(umap_df_filt, d_[, c("object_tag", umap_plot_markers)],
                         by = "object_tag", all.x = TRUE, sort = FALSE)

  cluster_plot <- ggplot2::ggplot(umap_df_filt,
                                   ggplot2::aes(UMAP1, UMAP2, color = factor(.data[[cluster_col]]))) +
    ggplot2::geom_point(size = 0.5, alpha = 0.9, shape = 16) +
    ggplot2::scale_color_manual(values = subcluster_colors) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "UMAP of nerve objects", color = "Subcluster") +
    cowplot::theme_cowplot() +
    ggplot2::theme(aspect.ratio = 1, axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))

  marker_plots <- lapply(umap_plot_markers, function(m) {
    ggplot2::ggplot(umap_df_filt, ggplot2::aes(UMAP1, UMAP2, color = .data[[m]])) +
      ggplot2::geom_point(size = 0.2, alpha = 0.8, shape = 16) +
      ggplot2::scale_color_gradientn(name = "Scaled\nMPI", limits = c(0, 1),
                                      values = scales::rescale(c(0.15, 1)),
                                      colours = cetcolor::cet_pal(9, name = "l8"),
                                      oob = scales::squish) +
      ggplot2::coord_equal() +
      ggplot2::labs(title = m) +
      cowplot::theme_cowplot() +
      ggplot2::theme(aspect.ratio = 1, axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
  })

  list(
    embedding_population = d,
    umap_df = umap_df,
    umap_df_filt = umap_df_filt,
    cluster_plot = cluster_plot,
    marker_plots = marker_plots,
    combined_plot = (patchwork::wrap_plots(cluster_plot,
                       patchwork::wrap_plots(marker_plots, ncol = 2, guides = "collect"),
                       ncol = 2) & ggplot2::theme(legend.position = "right"))
  )
}

#' WCSS (elbow) diagnostic for choosing k
#'
#' Loops `stats::kmeans(centers = k)$tot.withinss` over `kmin:kmax` and
#' converts to a *relative* drop (`-diff(tot_withinss) / tot_withinss`),
#' as shown in the "Relative WSS improvement per added cluster" plot.
#' `rel_drop` is `NA` for the first `k` by construction.
#'
#' @param df Data frame containing `cluster_cols`.
#' @param cluster_cols Character vector of marker columns to cluster on.
#' @param kmin,kmax Range of k to evaluate (default 2:12).
#' @param scale_for_kmeans Whether to `scale()` the marker matrix before
#'   `kmeans()` (default `TRUE`).
#' @param nstart,iter_max `stats::kmeans()` parameters (default 50, 100).
#'   Stage 1 clustering uses the `stats::kmeans()` default `nstart = 1`;
#'   this diagnostic uses a more thorough search to estimate a stable WSS
#'   curve.
#' @param seed Base seed; each k uses `seed + k` (default 123).
#' @param max_n If set and `nrow(df) > max_n`, downsample to `max_n` rows
#'   (seeded on `seed`) before clustering, for runtime on very large
#'   populations. `NULL` (default) runs on the full population.
#'
#' @return A tibble with `k`, `tot_withinss`, `delta`, `rel_drop`.
compute_elbow_wss <- function(df, cluster_cols, kmin = 2, kmax = 12,
                               scale_for_kmeans = TRUE, nstart = 50,
                               iter_max = 100, seed = 123, max_n = NULL) {
  stopifnot(all(cluster_cols %in% names(df)))

  x <- df %>%
    dplyr::select(dplyr::all_of(cluster_cols)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    tidyr::drop_na()

  if (!is.null(max_n) && nrow(x) > max_n) {
    set.seed(seed)
    x <- x %>% dplyr::slice_sample(n = max_n)
  }

  x_mat <- if (scale_for_kmeans) scale(as.matrix(x)) else as.matrix(x)
  ks <- kmin:kmax

  wss <- sapply(ks, function(k) {
    set.seed(seed + k)
    stats::kmeans(x_mat, centers = k, nstart = nstart, iter.max = iter_max)$tot.withinss
  })

  tibble::tibble(k = ks, tot_withinss = wss) %>%
    dplyr::mutate(delta = c(NA, diff(tot_withinss)), rel_drop = -delta / tot_withinss)
}

#' Silhouette diagnostic for choosing k
#'
#' Loops `stats::kmeans(centers = k)$cluster` over `kmin:kmax` and
#' computes the mean silhouette width (`cluster::silhouette()`) against
#' the same scaled marker matrix, as shown in the "Average silhouette
#' width" plot.
#'
#' @param df Data frame containing `cluster_cols`.
#' @param cluster_cols Character vector of marker columns to cluster on.
#' @param kmin,kmax Range of k to evaluate (default 2:12).
#' @param scale_for_kmeans Whether to `scale()` the marker matrix before
#'   `kmeans()` (default `TRUE`).
#' @param nstart `stats::kmeans()` parameter (default 25).
#' @param seed Base seed for downsampling and each k's `kmeans()` call
#'   (`seed + k`) (default 222, matching [cfg$clustering$seed]).
#' @param max_n If `nrow(df) > max_n`, downsample to `max_n` rows (seeded
#'   on `seed`) before clustering -- silhouette's O(n^2) distance matrix
#'   makes this far more necessary than for [compute_elbow_wss()]. Default
#'   2222; [plot_wcss_silhouette_diagnostics()] uses 8000.
#'
#' @return A tibble with `k`, `silhouette`.
calc_silhouette <- function(df, cluster_cols, kmin = 2, kmax = 12,
                             scale_for_kmeans = TRUE, nstart = 25,
                             seed = cfg$clustering$seed, max_n = 2222) {
  set.seed(seed)

  x <- df %>%
    dplyr::select(dplyr::all_of(cluster_cols)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    tidyr::drop_na()

  if (nrow(x) > max_n) x <- x %>% dplyr::slice_sample(n = max_n)

  x_mat <- if (scale_for_kmeans) scale(as.matrix(x)) else as.matrix(x)
  d <- stats::dist(x_mat)
  ks <- kmin:kmax

  sil <- sapply(ks, function(k) {
    set.seed(seed + k)
    cl <- stats::kmeans(x_mat, centers = k, nstart = nstart)$cluster
    mean(cluster::silhouette(cl, d)[, "sil_width"])
  })

  tibble::tibble(k = ks, silhouette = sil)
}

#' WCSS + silhouette diagnostic plots for one population
#'
#' Combines [compute_elbow_wss()] and [calc_silhouette()] into the
#' side-by-side "k selection diagnostics" figure. Run this ONCE on the
#' Stage 1 all-objects population (`step1$objs4pca`, i.e. the output of
#' [build_objs4pca()], to choose `cfg$clustering$k_og`) and ONCE on the Stage 2
#' ALL-nerve-candidate-objects population (`sub$nerve_objs`, i.e. the
#' full output of [subcluster_nerve_objects()] restricted to k=3 or
#' whichever single k you're diagnosing against -- NOT the smaller,
#' already-subcluster-3-excluded [run_nerve_object_kmeans_pipeline()]
#' `$nerve_obj_table`, to choose `cfg$clustering$k_sub_max`). Both are the SAME
#' populations `run_og_kmeans()`/`subcluster_nerve_objects()` themselves
#' cluster on -- this is a pre-clustering diagnostic, not a downstream one.
#'
#' @param df Data frame containing `quant_cols`.
#' @param quant_cols Character vector of marker columns. Defaults to [cfg$clustering$quant_cols].
#' @param kmin,kmax Range of k to evaluate (default 2:12).
#' @param wss_seed,sil_seed Base seeds for [compute_elbow_wss()]/
#'   [calc_silhouette()] (default 123/[cfg$clustering$seed]).
#' @param max_n Downsampling cap applied to BOTH diagnostics (default
#'   8000; `NULL` to run [compute_elbow_wss()] on the full population while
#'   still capping [calc_silhouette()] at its own default).
#' @param subtitle Plot subtitle (e.g. `"Unfiltered all-objects (Stroma)"`
#'   or `"Nerve-candidate objects"`).
#'
#' @return A `patchwork` object: WSS plot | silhouette plot.
plot_wcss_silhouette_diagnostics <- function(df, quant_cols = cfg$clustering$quant_cols,
                                              kmin = 2, kmax = 12,
                                              wss_seed = 123, sil_seed = cfg$clustering$seed,
                                              max_n = 8000, subtitle = NULL) {
  marker_cols <- setNames(quant_cols, vapply(quant_cols, .marker_shortname, character(1)))
  d <- df
  for (nm in names(marker_cols)) d[[nm]] <- d[[marker_cols[[nm]]]]

  elbow <- compute_elbow_wss(d, names(marker_cols), kmin = kmin, kmax = kmax,
                              seed = wss_seed, max_n = max_n)
  sil <- calc_silhouette(d, names(marker_cols), kmin = kmin, kmax = kmax,
                            seed = sil_seed, max_n = if (is.null(max_n)) 2222 else max_n)

  clean_theme <- ggplot2::theme_bw(base_size = 13) + ggplot2::theme(panel.grid = ggplot2::element_blank())

  p_wss <- ggplot2::ggplot(elbow, ggplot2::aes(k, rel_drop)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 2) + clean_theme +
    ggplot2::labs(title = "Relative WSS improvement per added cluster", subtitle = subtitle,
                  x = "k", y = "Relative drop in WSS") +
    ggplot2::scale_x_continuous(limits = c(0, 13), n.breaks = 7, expand = FALSE)

  p_sil <- ggplot2::ggplot(sil, ggplot2::aes(k, silhouette)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 2) + clean_theme +
    ggplot2::labs(title = "Average silhouette width", subtitle = subtitle,
                  x = "k", y = "Silhouette score") +
    ggplot2::scale_x_continuous(limits = c(0, 13), n.breaks = 7, expand = FALSE)

  p_wss + p_sil
}
