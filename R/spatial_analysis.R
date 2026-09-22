#' ============================================================
# --- Spatial analysis ---------------------------------------
#' ============================================================


#' Run spatial STaN enrichment analysis with Tweedie GLMMs
#'
#' Prepares image-level clinical/QC data, applies spatial band eligibility
#' criteria, fits one Tweedie GLMM per STaN/vessel subtype, performs planned
#' contrasts, applies BH correction, runs DHARMa diagnostics, and exports the
#' complete long-form statistics table.
#'
#' Reviewer-adjustable analysis parameters are stored in cfg$spatial_tweedie.
#' Graphical settings are intentionally handled elsewhere.
#'
#' @param patient_df Patient-level STaN summary containing
#'   `consolidated_im3s`, StudyID, clinical outcomes, and nerve_group.
#' @param image_df QC table containing image-level inclusion decisions.
#' @param cfg Project configuration list.
#' @param outcome Outcome to model: "BCG_failure" or "Progression".
#' @param band_file Path to nerve_vessel_area_by_distance_to_tumor_stroma_border.csv. Defaults to project path.
#' @param save_results Write long-form model statistics to TSV.
#' @param band_data Optional precomputed band-area table; replaces CSV reading only.
#'
#' @return List containing prepared data, fitted models, diagnostics, contrasts,
#'   and final long- and wide-form statistics tables.
#'
#' USAGE:
#'   res <- run_spatial_tweedie(
#'     patient_df,
#'     image_df,
#'     cfg,
#'     outcome = "Progression"
#'   )
run_spatial_tweedie <- function(patient_df, image_df, cfg,
                                outcome = "Progression",
                                band_file = NULL,
                                save_results = TRUE, band_data = NULL) {

  s <- cfg$spatial_tweedie
  patient_col <- "StudyID"; image_col <- "image_id"

  # ---- Validate reviewer-adjustable settings -------------------------------
  if (!outcome %in% names(s$outcome_levels))
    stop("No outcome levels configured for `", outcome, "`.")
  if (!s$p_adjust_scope %in% c("global", "subcluster"))
    stop("`p_adjust_scope` must be 'global' or 'subcluster'.")
  if (!s$weight_mode %in% c("none", "stroma", "nerve", "both"))
    stop("`weight_mode` must be 'none', 'stroma', 'nerve', or 'both'.")
  if (!s$weight_transform %in% c("sqrt", "linear"))
    stop("`weight_transform` must be 'sqrt' or 'linear'.")
  if (!is.null(s$nerve_weight_cap_quantile) &&
      (!is.numeric(s$nerve_weight_cap_quantile) ||
       length(s$nerve_weight_cap_quantile) != 1 ||
       s$nerve_weight_cap_quantile <= 0 ||
       s$nerve_weight_cap_quantile > 1))
    stop("`nerve_weight_cap_quantile` must be NULL or in (0, 1].")

  oc_levels <- s$outcome_levels[[outcome]]
  if (is.null(band_file))
    band_file <- file.path(cfg$processed_dir, "nerve_vessel_area_by_distance_to_tumor_stroma_border.csv")


  # ==========================================================================
  # 1. PREPARE IMAGE-LEVEL CLINICAL / QC DATA
  # ==========================================================================

  # Expand comma-separated image IDs so each image occupies one row.
  spatial_clin <- patient_df %>%
    tidyr::separate_rows(consolidated_im3s, sep = ", ", convert = TRUE) %>%
    dplyr::mutate(
      `Sample Name` = consolidated_im3s,
      image_id = consolidated_im3s,
      image_id = gsub(" |\\-|,", "_", image_id),
      image_id = gsub(".im3|\\]|\\[", "", image_id),
      image_id = gsub("\\#", "slide", image_id)
    )

  # Published inclusion flags are row-level decisions, not a scalar test of
  # the entire column. Both archived image IDs and source filenames are supported.
  keep <- image_df %>%
    dplyr::filter(include_image_final_spatial_cohort %in% TRUE)
  if ("im3_filename" %in% names(keep)) {
    spatial_clin <- spatial_clin %>%
      dplyr::filter(`Sample Name` %in% keep$im3_filename)
  } else {
    spatial_clin <- spatial_clin %>%
      dplyr::filter(image_id %in% keep$image_id)
  }

  # Spatial outcome comparisons use the study's adequate-BCG analytical cohort.
  # Require documented Yes (missing adequacy does not establish eligibility).
  # Apply this to analysis data, not to upstream spatial measurements.
  spatial_clin <- spatial_clin %>%
    dplyr::filter(Adequate_BCG == "Yes")


  # ==========================================================================
  # 2. LOAD + PREPARE BORDER-BAND DATA
  # ==========================================================================

  if (is.null(band_data)) {
    band_data <- readr::read_delim(
      band_file, delim = ",", trim_ws = TRUE,
      col_names = TRUE, show_col_types = FALSE
    )
  }
  nerve_border_bands <- as.data.frame(band_data) %>%

    # Only retain bands within the configured spatial radius.
    dplyr::filter(band_outer_um <= s$max_radius_um) %>%

    # Assign nerve/vessel subtype identifiers and harmonize band labels.
    dplyr::mutate(
      subcluster = ifelse(
        layer_name == "nerve",
        paste0("subcluster_", subtype_code),
        "Vessel"
      ),
      subtype_label = dplyr::case_when(
        subcluster == "subcluster_1" ~ "SP+/VGLUT1+ sc1",
        subcluster == "subcluster_2" ~ "TH.mid sc2",
        subcluster == "subcluster_4" ~ "VACHT+ sc4",
        subcluster == "subcluster_5" ~ "TH.high sc5",
        subcluster == "subcluster_6" ~ "SYP.high sc6",
        subcluster == "subcluster_7" ~ "TH-/VACHT-/SP-/VGLUT1- sc7",
        subcluster == "Vessel" ~ "Vessel",
        TRUE ~ NA_character_
      ),
      band_label = dplyr::if_else(
        band_label == "200-300 um", "201-300 um", band_label
      )
    ) %>%

    # Calculate radius-level quantities separately for each image/subtype.
    dplyr::group_by(image_id, subtype_code, subcluster, subtype_label) %>%
    dplyr::mutate(
      band_label = factor(band_label, levels = s$band_levels),
      max_radius = max(band_outer_um),

      # Stroma in the outermost available band determines whether the image
      # has sufficient tissue coverage out to the requested radius.
      stroma_area_um2_at_max_radius_band =
        stroma_area_um2[match(max_radius, band_outer_um)],

      radius_area_um2 = sum(stroma_area_um2),

      # Store actual nerve/vessel area rather than reconstructing it later
      # from a percentage-valued radius_density.
      radius_nerve_area_um2 = sum(area_um2),

      radius_density = (radius_nerve_area_um2 / radius_area_um2) * 100,
      band_density = (area_um2 / stroma_area_um2) * 100,

      # Radius-normalized enrichment: local band density / overall density
      # within the complete spatial radius.
      band_enrichment = dplyr::if_else(
        radius_density > 0,
        band_density / radius_density,
        NA_real_
      ),
      log2_band_enrichment = log2(band_enrichment + 0.001),
      smallest_band_area = min(stroma_area_um2)
    ) %>%

    # Apply spatial coverage eligibility ONCE here.
    # These filters therefore do not need to be repeated before model fitting.
    dplyr::filter(
      max_radius == s$max_radius_um,
      stroma_area_um2_at_max_radius_band >= s$min_outer_band_area_um2,
      smallest_band_area >= s$min_band_area_um2
    ) %>%
    dplyr::ungroup() %>%

    # Attach patient-level clinical variables.
    dplyr::left_join(
      spatial_clin %>%
        dplyr::select(
          image_id, StudyID, `Sample Name`,
          BCG_failure, Progression, nerve_group
        ),
      by = "image_id"
    )

  # Optional restriction to a patient-level STaN abundance group.
  # Current manuscript analysis uses STaN-high tumors only.
  if (!is.null(s$nerve_group_keep))
    nerve_border_bands <- nerve_border_bands %>%
    dplyr::filter(nerve_group %in% s$nerve_group_keep)

  if (!nrow(nerve_border_bands))
    stop("No border-band observations remain after spatial eligibility filtering.")


  # Image list retained for downstream spatial visualizations if needed.
  bands_keep <- nerve_border_bands %>%
    dplyr::group_by(
      image_id, `Sample Name`, BCG_failure,
      Progression, nerve_group
    ) %>%
    dplyr::summarise(
      stroma_area_um2 = sum(stroma_area_um2),
      StudyID = dplyr::first(StudyID),
      .groups = "drop"
    )


  # ==========================================================================
  # 3. PREPARE TWEE­DIE MODEL DATA
  # ==========================================================================

  required <- c(
    patient_col, image_col, outcome,
    "subcluster", "subtype_label", "band_label",
    "area_um2", "stroma_area_um2", "radius_area_um2",
    "radius_nerve_area_um2", "band_density", "radius_density"
  )
  missing <- setdiff(required, names(nerve_border_bands))
  if (length(missing))
    stop("Missing required columns: ", paste(missing, collapse = ", "))

  model_df <- nerve_border_bands %>%
    dplyr::filter(
      !is.na(.data[[patient_col]]),
      !is.na(.data[[image_col]]),
      !is.na(.data[[outcome]]),
      !is.na(subcluster),
      !is.na(subtype_label),
      !is.na(band_label),

      # Band-level zeros are valid Tweedie observations.
      is.finite(band_density), band_density >= 0,

      # An overall radius density >0 is required for the offset/enrichment
      # formulation and for nerve-area weighting.
      is.finite(radius_density), radius_density > 0,
      is.finite(radius_nerve_area_um2),
      radius_nerve_area_um2 > s$min_radius_nerve_area_um2
    )

  # Optional subtype exclusion, currently excluding cluster 1.
  if (!is.null(s$exclude_subclusters))
    model_df <- model_df %>%
    dplyr::filter(!subcluster %in% s$exclude_subclusters)

  model_df <- model_df %>%
    dplyr::mutate(
      patient_model = factor(.data[[patient_col]]),
      image_model = factor(.data[[image_col]]),
      outcome_model = factor(.data[[outcome]], levels = oc_levels),
      subcluster = factor(subcluster),
      subtype_label = factor(subtype_label),
      band_label = factor(band_label, levels = s$band_levels),
      log_radius_density = log(radius_density)
    ) %>%
    droplevels()

  if (!nrow(model_df))
    stop("No observations remain after model filtering.")
  if (nlevels(model_df$outcome_model) != 2)
    stop(
      "Outcome must contain both configured levels after filtering: ",
      paste(oc_levels, collapse = ", ")
    )


  # Lookup used to restore human-readable subtype labels to result tables.
  subtype_lookup <- model_df %>%
    dplyr::distinct(subcluster, subtype_label) %>%
    dplyr::transmute(
      subcluster = as.character(subcluster),
      subtype_label = as.character(subtype_label)
    )

  if (anyDuplicated(subtype_lookup$subcluster))
    stop("At least one subcluster maps to multiple subtype labels.")

  add_subtype <- function(x, sc) {
    lab <- subtype_lookup$subtype_label[
      match(sc, subtype_lookup$subcluster)
    ]
    dplyr::mutate(
      x,
      subcluster = sc,
      subtype_label = lab,
      .before = 1
    )
  }


  # ==========================================================================
  # 4. FIT ONE TWEE­DIE GLMM PER SUBTYPE
  # ==========================================================================

  # band_density / radius_density represents radius-normalized enrichment.
  # log(radius_density) is therefore included as an offset.
  #
  # Images are nested within patients to account for multiple images/patient.
  model_formula <- band_density ~ outcome_model * band_label +
    offset(log_radius_density) +
    (1 | patient_model / image_model)

  fit_one <- function(dat) {
    dat <- droplevels(dat)

    # Do not attempt models without enough levels/clusters for estimation.
    if (nlevels(dat$outcome_model) != 2 ||
        nlevels(dat$band_label) < 2 ||
        dplyr::n_distinct(dat$patient_model) < 2 ||
        dplyr::n_distinct(dat$image_model) < 2)
      return(NULL)

    nerve_area <- dat$radius_nerve_area_um2

    # Prevent the largest nerve-rich observations from dominating weights.
    if (!is.null(s$nerve_weight_cap_quantile)) {
      cap <- stats::quantile(
        nerve_area,
        s$nerve_weight_cap_quantile,
        na.rm = TRUE, names = FALSE, type = 8
      )
      nerve_area <- pmin(nerve_area, cap)
    }

    wt <- function(x)
      if (s$weight_transform == "sqrt") sqrt(x) else x

    raw_weight <- switch(
      s$weight_mode,
      none   = rep(1, nrow(dat)),
      stroma = wt(dat$stroma_area_um2),
      nerve  = wt(nerve_area),
      both   = wt(dat$stroma_area_um2) * wt(nerve_area)
    )

    # Scaling weights to mean 1 retains relative weighting while making
    # magnitude comparable across subtype-specific models.
    dat$analysis_weight <- raw_weight / mean(raw_weight, na.rm = TRUE)

    if (any(!is.finite(dat$analysis_weight)) ||
        any(dat$analysis_weight <= 0))
      stop(
        "Invalid weights for ",
        paste(unique(as.character(dat$subcluster)), collapse = ", ")
      )

    tryCatch(
      glmmTMB::glmmTMB(
        model_formula,
        family = glmmTMB::tweedie(link = "log"),
        weights = analysis_weight,
        data = dat
      ),
      error = \(e) {
        warning(
          "Model failed for ",
          unique(as.character(dat$subcluster)),
          ": ", e$message
        )
        NULL
      }
    )
  }

  subsets <- split(model_df, model_df$subcluster, drop = TRUE)

  model_list <- purrr::imap(subsets, \(dat, sc) {
    message("Fitting ", sc, " using ", s$weight_mode, " weights...")
    fit_one(dat)
  }) %>%
    purrr::compact()

  if (!length(model_list))
    stop("No subtype models successfully fitted.")


  # ==========================================================================
  # 5. CONTRAST HELPERS
  # ==========================================================================

  # For each band: compare that band with the equally weighted mean of all
  # remaining distance bands.
  one_vs_rest <- function(bands) {
    out <- lapply(seq_along(bands), \(i) {
      x <- rep(-1 / (length(bands) - 1), length(bands))
      x[i] <- 1
      x
    })
    names(out) <- paste0(bands, "_vs_other_bands")
    out
  }

  # Tweedie uses a natural-log link, so `estimate` is a log ratio and
  # exp(estimate) is the corresponding enrichment ratio.
  add_ratio_columns <- function(x) {
    lcl <- intersect(c("lower.CL", "asymp.LCL"), names(x))[1]
    ucl <- intersect(c("upper.CL", "asymp.UCL"), names(x))[1]

    dplyr::mutate(
      x,
      log_ratio = estimate,
      ratio = exp(estimate),
      ratio_lower = if (length(lcl)) exp(.data[[lcl]]) else NA_real_,
      ratio_upper = if (length(ucl)) exp(.data[[ucl]]) else NA_real_
    )
  }


  # ==========================================================================
  # 6. MODEL / CONVERGENCE SUMMARY
  # ==========================================================================

  model_diagnostics_df <- purrr::imap_dfr(model_list, \(m, sc) {
    mf <- stats::model.frame(m)

    add_subtype(
      tibble::tibble(
        n_observations = stats::nobs(m),
        n_patients = dplyr::n_distinct(mf$patient_model),
        n_images = dplyr::n_distinct(mf$image_model),
        convergence_code = m$fit$convergence,
        positive_definite_hessian = isTRUE(m$sdr$pdHess),
        tweedie_power = unname(
          stats::plogis(m$fit$parfull["psi"]) + 1
        ),
        weight_mode = s$weight_mode,
        weight_transform = s$weight_transform,
        nerve_weight_cap_quantile =
          if (is.null(s$nerve_weight_cap_quantile))
            NA_real_ else s$nerve_weight_cap_quantile
      ),
      sc
    )
  })


  # ==========================================================================
  # 7. OMNIBUS TESTS + ESTIMATED ENRICHMENT
  # ==========================================================================

  omnibus_tests_df <- purrr::imap_dfr(
    model_list,
    \(m, sc)
    as.data.frame(emmeans::joint_tests(m)) %>%
      add_subtype(sc)
  )

  # Setting the offset to zero returns the modeled band/radius enrichment
  # rather than predictions at a particular observed radius density.
  estimated_means_df <- purrr::imap_dfr(model_list, \(m, sc) {
    z <- as.data.frame(
      emmeans::emmeans(
        m,
        ~ outcome_model * band_label,
        offset = 0,
        type = "response"
      )
    )

    est <- intersect(c("response", "rate", "prob", "mu", "emmean"), names(z))[1]
    lcl <- intersect(c("asymp.LCL", "lower.CL"), names(z))[1]
    ucl <- intersect(c("asymp.UCL", "upper.CL"), names(z))[1]

    if (!length(est) || !length(lcl) || !length(ucl))
      stop("Could not identify emmeans response columns for ", sc)

    z %>%
      dplyr::mutate(
        estimated_enrichment = .data[[est]],
        enrichment_lower = .data[[lcl]],
        enrichment_upper = .data[[ucl]]
      ) %>%
      add_subtype(sc)
  })


  # ==========================================================================
  # 8. PLANNED CONTRASTS
  # ==========================================================================

  # Outcome comparison within each distance band.
  stats_between_outcomes_df <- purrr::imap_dfr(model_list, \(m, sc) {
    emmeans::contrast(
      emmeans::emmeans(
        m, ~ outcome_model | band_label, offset = 0
      ),
      method = "revpairwise",
      adjust = "none"
    ) %>%
      summary(infer = c(TRUE, TRUE), adjust = "none") %>%
      as.data.frame() %>%
      add_ratio_columns() %>%
      add_subtype(sc)
  })

  # Each band versus the mean of all other bands, separately by outcome.
  stats_within_outcomes_df <- purrr::imap_dfr(model_list, \(m, sc) {
    bands <- levels(stats::model.frame(m)$band_label)

    emmeans::contrast(
      emmeans::emmeans(
        m, ~ band_label | outcome_model, offset = 0
      ),
      method = one_vs_rest(bands),
      adjust = "none"
    ) %>%
      summary(infer = c(TRUE, TRUE), adjust = "none") %>%
      as.data.frame() %>%
      add_ratio_columns() %>%
      dplyr::mutate(
        test_band = sub("_vs_other_bands$", "", contrast)
      ) %>%
      add_subtype(sc)
  })

  # Compare each band-vs-rest contrast between the two outcomes:
  # i.e. does the spatial enrichment pattern differ by clinical outcome?
  stats_spatial_pattern_between_outcomes_df <-
    purrr::imap_dfr(model_list, \(m, sc) {

      bands <- levels(stats::model.frame(m)$band_label)

      band_contrasts <- emmeans::contrast(
        emmeans::emmeans(
          m, ~ band_label | outcome_model, offset = 0
        ),
        method = one_vs_rest(bands),
        by = "outcome_model",
        adjust = "none"
      )

      z <- emmeans::contrast(
        band_contrasts,
        method = "revpairwise",
        by = "contrast",
        adjust = "none"
      ) %>%
        summary(infer = c(TRUE, TRUE), adjust = "none") %>%
        as.data.frame()

      outcome_col <- intersect(
        c("contrast1", "outcome_model"), names(z)
      )[1]

      if (!length(outcome_col))
        stop("Could not identify outcome contrast column for ", sc)

      z %>%
        add_ratio_columns() %>%
        dplyr::mutate(
          test_band = sub("_vs_other_bands$", "", contrast),
          outcome_contrast = .data[[outcome_col]]
        ) %>%
        add_subtype(sc)
    })


  # ==========================================================================
  # 9. BH MULTIPLE-TESTING CORRECTION
  # ==========================================================================

  if (s$p_adjust_scope == "subcluster") {

    stats_between_outcomes_df <- stats_between_outcomes_df %>%
      dplyr::group_by(subcluster, subtype_label) %>%
      dplyr::mutate(
        padj_between_outcomes = p.adjust(p.value, "BH")
      ) %>%
      dplyr::ungroup()

    stats_within_outcomes_df <- stats_within_outcomes_df %>%
      dplyr::group_by(subcluster, subtype_label, outcome_model) %>%
      dplyr::mutate(
        padj_within_outcome = p.adjust(p.value, "BH")
      ) %>%
      dplyr::ungroup()

    stats_spatial_pattern_between_outcomes_df <-
      stats_spatial_pattern_between_outcomes_df %>%
      dplyr::group_by(subcluster, subtype_label) %>%
      dplyr::mutate(
        padj_spatial_pattern_between_outcomes = p.adjust(p.value, "BH")
      ) %>%
      dplyr::ungroup()

  } else {

    stats_between_outcomes_df <- stats_between_outcomes_df %>%
      dplyr::mutate(
        padj_between_outcomes = p.adjust(p.value, "BH")
      )

    stats_within_outcomes_df <- stats_within_outcomes_df %>%
      dplyr::mutate(
        padj_within_outcome = p.adjust(p.value, "BH")
      )

    stats_spatial_pattern_between_outcomes_df <-
      stats_spatial_pattern_between_outcomes_df %>%
      dplyr::mutate(
        padj_spatial_pattern_between_outcomes =
          p.adjust(p.value, "BH")
      )
  }


  # ==========================================================================
  # 10. SAMPLE COUNTS
  # ==========================================================================

  sample_counts_by_outcome_df <- model_df %>%
    dplyr::group_by(
      subcluster, subtype_label, outcome_model
    ) %>%
    dplyr::summarise(
      n_patients = dplyr::n_distinct(patient_model),
      n_images = dplyr::n_distinct(image_model),
      n_rows = dplyr::n(),
      n_zero_bands = sum(band_density == 0, na.rm = TRUE),
      prop_zero_bands = mean(band_density == 0, na.rm = TRUE),
      median_radius_nerve_area_um2 =
        median(radius_nerve_area_um2, na.rm = TRUE),
      .groups = "drop"
    )

  sample_counts_by_outcome_band_df <- model_df %>%
    dplyr::group_by(
      subcluster, subtype_label, outcome_model, band_label
    ) %>%
    dplyr::summarise(
      n_patients = dplyr::n_distinct(patient_model),
      n_images = dplyr::n_distinct(image_model),
      n_rows = dplyr::n(),
      n_zero_bands = sum(band_density == 0, na.rm = TRUE),
      prop_zero_bands = mean(band_density == 0, na.rm = TRUE),
      median_band_stroma_area_um2 =
        median(stroma_area_um2, na.rm = TRUE),
      total_band_stroma_area_um2 =
        sum(stroma_area_um2, na.rm = TRUE),
      median_radius_nerve_area_um2 =
        median(radius_nerve_area_um2, na.rm = TRUE),
      .groups = "drop"
    )


  # ==========================================================================
  # 11. BUILD FINAL WIDE + LONG RESULT TABLES
  # ==========================================================================

  sample_counts_wide_df <- sample_counts_by_outcome_band_df %>%
    dplyr::mutate(
      outcome_model = paste0("outcome_", outcome_model)
    ) %>%
    tidyr::pivot_wider(
      names_from = outcome_model,
      values_from = c(
        n_patients, n_images, n_rows,
        n_zero_bands, prop_zero_bands,
        median_band_stroma_area_um2,
        total_band_stroma_area_um2,
        median_radius_nerve_area_um2
      ),
      names_glue = "{.value}_{outcome_model}"
    )

  estimated_means_wide_df <- estimated_means_df %>%
    dplyr::mutate(
      outcome_model = paste0("outcome_", outcome_model)
    ) %>%
    dplyr::select(
      subcluster, subtype_label, band_label, outcome_model,
      estimated_enrichment, enrichment_lower,
      enrichment_upper, SE
    ) %>%
    tidyr::pivot_wider(
      names_from = outcome_model,
      values_from = c(
        estimated_enrichment,
        enrichment_lower,
        enrichment_upper,
        SE
      ),
      names_glue = "{.value}_{outcome_model}"
    )

  within_outcome_wide_df <- stats_within_outcomes_df %>%
    dplyr::mutate(
      outcome_model = paste0("outcome_", outcome_model)
    ) %>%
    dplyr::transmute(
      subcluster, subtype_label, test_band, outcome_model,
      band_vs_rest_log_ratio = log_ratio,
      band_vs_rest_ratio = ratio,
      band_vs_rest_ratio_lower = ratio_lower,
      band_vs_rest_ratio_upper = ratio_upper,
      pval_within_outcome = p.value,
      padj_within_outcome
    ) %>%
    tidyr::pivot_wider(
      names_from = outcome_model,
      values_from = c(
        band_vs_rest_log_ratio,
        band_vs_rest_ratio,
        band_vs_rest_ratio_lower,
        band_vs_rest_ratio_upper,
        pval_within_outcome,
        padj_within_outcome
      ),
      names_glue = "{.value}_{outcome_model}"
    )

  pvals_model_stats <- estimated_means_wide_df %>%
    dplyr::left_join(
      stats_between_outcomes_df %>%
        dplyr::transmute(
          subcluster, subtype_label, band_label,
          outcome_contrast = contrast,
          outcome_log_ratio = log_ratio,
          outcome_enrichment_ratio = ratio,
          outcome_ratio_lower = ratio_lower,
          outcome_ratio_upper = ratio_upper,
          pval_between_outcomes = p.value,
          padj_between_outcomes
        ),
      by = c("subcluster", "subtype_label", "band_label")
    ) %>%
    dplyr::left_join(
      stats_spatial_pattern_between_outcomes_df %>%
        dplyr::transmute(
          subcluster, subtype_label, test_band,
          spatial_pattern_outcome_contrast = outcome_contrast,
          spatial_pattern_log_ratio = log_ratio,
          spatial_pattern_ratio = ratio,
          spatial_pattern_ratio_lower = ratio_lower,
          spatial_pattern_ratio_upper = ratio_upper,
          pval_spatial_pattern_between_outcomes = p.value,
          padj_spatial_pattern_between_outcomes
        ),
      by = c(
        "subcluster", "subtype_label",
        "band_label" = "test_band"
      )
    ) %>%
    dplyr::left_join(
      within_outcome_wide_df,
      by = c(
        "subcluster", "subtype_label",
        "band_label" = "test_band"
      )
    ) %>%
    dplyr::left_join(
      sample_counts_wide_df,
      by = c("subcluster", "subtype_label", "band_label")
    ) %>%
    dplyr::arrange(subcluster, band_label)

  pvals_model_stats_long <- estimated_means_df %>%
    dplyr::transmute(
      subcluster, subtype_label, band_label,
      outcome = as.character(outcome_model),
      estimated_enrichment,
      enrichment_lower,
      enrichment_upper,
      enrichment_SE = SE
    ) %>%
    dplyr::left_join(
      stats_within_outcomes_df %>%
        dplyr::transmute(
          subcluster, subtype_label,
          band_label = test_band,
          outcome = as.character(outcome_model),
          band_vs_rest_log_ratio = log_ratio,
          band_vs_rest_ratio = ratio,
          band_vs_rest_ratio_lower = ratio_lower,
          band_vs_rest_ratio_upper = ratio_upper,
          pval_within_outcome = p.value,
          padj_within_outcome
        ),
      by = c(
        "subcluster", "subtype_label",
        "band_label", "outcome"
      )
    ) %>%
    dplyr::left_join(
      sample_counts_by_outcome_band_df %>%
        dplyr::transmute(
          subcluster, subtype_label, band_label,
          outcome = as.character(outcome_model),
          n_patients, n_images, n_rows,
          n_zero_bands, prop_zero_bands,
          median_band_stroma_area_um2,
          total_band_stroma_area_um2,
          median_radius_nerve_area_um2
        ),
      by = c(
        "subcluster", "subtype_label",
        "band_label", "outcome"
      )
    ) %>%
    dplyr::left_join(
      stats_between_outcomes_df %>%
        dplyr::transmute(
          subcluster, subtype_label, band_label,
          outcome_contrast = contrast,
          outcome_log_ratio = log_ratio,
          outcome_enrichment_ratio = ratio,
          outcome_ratio_lower = ratio_lower,
          outcome_ratio_upper = ratio_upper,
          pval_between_outcomes = p.value,
          padj_between_outcomes
        ),
      by = c("subcluster", "subtype_label", "band_label")
    ) %>%
    dplyr::left_join(
      stats_spatial_pattern_between_outcomes_df %>%
        dplyr::transmute(
          subcluster, subtype_label,
          band_label = test_band,
          spatial_pattern_outcome_contrast = outcome_contrast,
          spatial_pattern_log_ratio = log_ratio,
          spatial_pattern_ratio = ratio,
          spatial_pattern_ratio_lower = ratio_lower,
          spatial_pattern_ratio_upper = ratio_upper,
          pval_spatial_pattern_between_outcomes = p.value,
          padj_spatial_pattern_between_outcomes
        ),
      by = c("subcluster", "subtype_label", "band_label")
    ) %>%
    dplyr::arrange(
      subcluster,
      band_label,
      factor(outcome, levels = oc_levels)
    ) %>%
    dplyr::mutate(
      band_label = factor(band_label, levels = s$band_levels),
      outcome = factor(outcome, levels = oc_levels),
      weight_mode = s$weight_mode,
      weight_transform = s$weight_transform,
      nerve_weight_cap_quantile =
        if (is.null(s$nerve_weight_cap_quantile))
          NA_real_ else s$nerve_weight_cap_quantile
    )


  # ==========================================================================
  # 12. DHARMa RESIDUAL DIAGNOSTICS
  # ==========================================================================

  dharma_tests_df <- purrr::imap_dfr(model_list, \(m, sc) {
    sim <- DHARMa::simulateResiduals(
      m,
      n = s$dharma_n_sim,
      plot = FALSE,
      seed = s$dharma_seed
    )

    add_subtype(
      tibble::tibble(
        weight_mode = s$weight_mode,
        uniformity_p =
          DHARMa::testUniformity(sim, plot = FALSE)$p.value,
        dispersion_p =
          DHARMa::testDispersion(sim, plot = FALSE)$p.value,
        zero_inflation_p =
          DHARMa::testZeroInflation(sim, plot = FALSE)$p.value,
        outlier_p =
          DHARMa::testOutliers(sim, plot = FALSE)$p.value
      ),
      sc
    )
  })


  # ==========================================================================
  # 13. EXPORT
  # ==========================================================================

  out_tsv <- NULL

  if (save_results) {
    out_dir <- file.path(cfg$path, cfg$paths$results$spatial, "tables")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    group_tag <- if (is.null(s$nerve_group_keep)) {
      "allnervegroups"
    } else {
      paste0(
        "nervegroup",
        paste(s$nerve_group_keep, collapse = "-")
      )
    }

    out_tsv <- file.path(
      out_dir,
      paste0(
        "spatial_tweedie_",
        group_tag,
        "_maxrad", s$max_radius_um,
        "_", outcome, ".tsv"
      )
    )

    readr::write_delim(
      pvals_model_stats_long,
      file = out_tsv,
      delim = "\t"
    )

    message("Saved: ", normalizePath(out_tsv, mustWork = FALSE))
  }


  # Return everything needed for figures, sensitivity analyses, or inspection.
  invisible(list(
    spatial_clin = spatial_clin,
    nerve_border_bands = nerve_border_bands,
    bands_keep = bands_keep,
    model_data = model_df,
    models = model_list,

    estimated_means = estimated_means_df,
    stats_between_outcomes = stats_between_outcomes_df,
    stats_within_outcomes = stats_within_outcomes_df,
    stats_spatial_pattern_between_outcomes =
      stats_spatial_pattern_between_outcomes_df,

    omnibus_tests = omnibus_tests_df,
    model_diagnostics = model_diagnostics_df,
    dharma_tests = dharma_tests_df,

    sample_counts = sample_counts_by_outcome_df,
    sample_counts_by_band = sample_counts_by_outcome_band_df,

    stats_wide = pvals_model_stats,
    stats_long = pvals_model_stats_long,

    output_file = out_tsv
  ))
}

###' ------------------------------------------------
# helper to generate SPLIT BUBBLE PLOTS to visualize differential band enrichment
###' ------------------------------------------------

.half_disc_polygon <- function(cx, cy, r, side = c("left", "right"), n = 40) {
  side <- match.arg(side)
  if (!is.finite(cx) || !is.finite(cy) || !is.finite(r) || r <= 0) return(NULL)
  theta <- if (side == "left") seq(pi / 2, 3 * pi / 2, length.out = n) else seq(-pi / 2, pi / 2, length.out = n)
  data.frame(x = c(cx, cx + r * cos(theta), cx), y = c(cy, cy + r * sin(theta), cy))
}

.full_circle_polygon <- function(cx, cy, r, n = 60) {
  if (!is.finite(cx) || !is.finite(cy) || !is.finite(r) || r <= 0) return(NULL)
  theta <- seq(0, 2 * pi, length.out = n + 1)
  data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

.build_split_circle_polygons <- function(d, x_col, y_col, id_cols, outcome_col,
                                         outcome_levels, z_col, magnitude_col,
                                         max_radius = 0.42, n_arc = 40) {
  d <- as.data.frame(d)
  x_levels <- unique(as.character(d[[x_col]]))
  y_levels <- unique(as.character(d[[y_col]]))
  d$cx_numeric <- match(as.character(d[[x_col]]), x_levels)
  d$cy_numeric <- match(as.character(d[[y_col]]), y_levels)

  cell_keys <- paste(d$cx_numeric, d$cy_numeric, sep = "_")
  max_mag_per_cell <- tapply(abs(d[[magnitude_col]]), cell_keys, max, na.rm = TRUE)
  d$shared_mag <- as.numeric(max_mag_per_cell[cell_keys])

  global_max <- max(d$shared_mag, na.rm = TRUE)
  if (!is.finite(global_max) || global_max <= 0) return(tibble::tibble())
  .scale_r <- function(v) max_radius * sqrt(v / global_max)

  d$outcome_str <- as.character(d[[outcome_col]])
  rows <- vector("list", nrow(d))

  for (i in seq_len(nrow(d))) {
    row <- d[i, , drop = FALSE]
    if (!row$outcome_str %in% outcome_levels) next
    side <- if (row$outcome_str == as.character(outcome_levels[1])) "left" else "right"

    poly <- .half_disc_polygon(
      cx = row$cx_numeric,
      cy = row$cy_numeric,
      r = .scale_r(row$shared_mag),
      side = side,
      n = n_arc
    )
    if (is.null(poly)) next

    poly$poly_id <- paste(c(as.character(unlist(row[id_cols])), row$outcome_str), collapse = "_")
    poly$z_value <- row[[z_col]]
    poly[[x_col]] <- row[[x_col]]
    poly[[y_col]] <- row[[y_col]]
    rows[[i]] <- poly
  }

  dplyr::bind_rows(rows)
}

# build rings
.build_significance_rings <- function(d, x_col, y_col, id_cols, magnitude_col,
                                      difference_p_col, p_threshold = 0.05,
                                      max_radius = 0.42, n_arc = 60) {
  d <- as.data.frame(d)
  if (!difference_p_col %in% names(d)) stop("Column '", difference_p_col, "' was not found.")

  x_levels <- unique(as.character(d[[x_col]]))
  y_levels <- unique(as.character(d[[y_col]]))

  sig_cells <- d %>%
    dplyr::filter(!is.na(.data[[difference_p_col]]), .data[[difference_p_col]] < p_threshold) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(id_cols, x_col, y_col))))

  if (!nrow(sig_cells)) return(tibble::tibble())

  sig_cells$cx_numeric <- match(as.character(sig_cells[[x_col]]), x_levels)
  sig_cells$cy_numeric <- match(as.character(sig_cells[[y_col]]), y_levels)

  all_x <- match(as.character(d[[x_col]]), x_levels)
  all_y <- match(as.character(d[[y_col]]), y_levels)
  cell_keys_all <- paste(all_x, all_y, sep = "_")
  max_mag_per_cell <- tapply(abs(d[[magnitude_col]]), cell_keys_all, max, na.rm = TRUE)

  global_max <- max(max_mag_per_cell, na.rm = TRUE)
  if (!is.finite(global_max) || global_max <= 0) return(tibble())
  .scale_r <- function(v) max_radius * sqrt(v / global_max)

  rows <- vector("list", nrow(sig_cells))

  for (i in seq_len(nrow(sig_cells))) {
    row <- sig_cells[i, , drop = FALSE]
    cell_key <- paste(row$cx_numeric, row$cy_numeric, sep = "_")
    cell_mag <- as.numeric(max_mag_per_cell[cell_key])
    if (!is.finite(cell_mag)) next

    poly <- .full_circle_polygon(
      cx = row$cx_numeric,
      cy = row$cy_numeric,
      r = .scale_r(cell_mag),
      n = n_arc
    )
    if (is.null(poly)) next

    poly$poly_id <- paste(as.character(unlist(row[id_cols])), collapse = "_")
    poly[[x_col]] <- row[[x_col]]
    poly[[y_col]] <- row[[y_col]]
    rows[[i]] <- poly
  }

  dplyr::bind_rows(rows)
}

# custom legends for half moons and sig rings
.make_half_circle_legend <- function(outcome_levels = c("No", "Yes"),
                                     title = "Outcome",
                                     sig_label = "FDR < 0.05",
                                     base_size = 13) {
  yy <- c(.40, .23, .06)
  labs <- tibble::tibble(x = .105, y = yy, label = c(outcome_levels, sig_label))
  polys <- dplyr::bind_rows(
    dplyr::mutate(.half_disc_polygon(.045, yy[1], .032, "left", 50), id = 1, type = "half"),
    dplyr::mutate(.half_disc_polygon(.045, yy[2], .032, "right", 50), id = 2, type = "half"),
    dplyr::mutate(.full_circle_polygon(.045, yy[3], .032, 70), id = 3, type = "ring")
  )

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = dplyr::filter(polys, type == "half"),
      ggplot2::aes(x = x, y = y, group = id),
      fill = "grey65", colour = NA
    ) +
    ggplot2::geom_polygon(
      data = dplyr::filter(polys, type == "ring"),
      ggplot2::aes(x = x, y = y, group = id),
      fill = NA, colour = "black", linewidth = .5
    ) +
    ggplot2::geom_text(
      data = labs,
      ggplot2::aes(x = x, y = y, label = label),
      hjust = 0, size = base_size / ggplot2::.pt
    ) +
    ggplot2::annotate(
      "text", x = 0, y = .565, label = title,
      hjust = 0, fontface = "plain", size = base_size / ggplot2::.pt
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, .62), ylim = c(0, .59),
      expand = FALSE, clip = "off"
    ) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0), aspect.ratio = 1)
}


# generate bubble plot!
plot_split_circle_grid <- function(pvals_df, x_col, y_col, outcome_col,
                                   outcome_levels = c("No", "Yes"),
                                   z_col = "log2FC",
                                   magnitude_col = "log2FC",
                                   difference_p_col = "padj_difference",
                                   difference_p_threshold = 0.05,
                                   title = NULL,
                                   x_lab = "Distance",
                                   y_lab = "Subcluster",
                                   fill_lab = NULL,
                                   split_lab = NULL,
                                   magnitude_lab = NULL,
                                   low = "blue2",
                                   mid = "white",
                                   high = "red2",
                                   max_radius = 0.42) {
  d <- as.data.frame(pvals_df)

  required_cols <- c(x_col, y_col, outcome_col, z_col, magnitude_col, difference_p_col)
  missing_cols <- dplyr::setdiff(required_cols, names(d))
  if (length(missing_cols)) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

  d[[outcome_col]] <- as.character(d[[outcome_col]])
  d <- d[d[[outcome_col]] %in% outcome_levels, , drop = FALSE]

  d_valid <- d[
    !is.na(d[[z_col]]) &
      !is.na(d[[magnitude_col]]) &
      is.finite(d[[z_col]]) &
      is.finite(d[[magnitude_col]]),
    ,
    drop = FALSE
  ]

  if (!nrow(d_valid)) return(ggplot2::ggplot() + cowplot::theme_cowplot())

  col_order <- if (is.factor(pvals_df[[x_col]])) levels(droplevels(pvals_df[[x_col]])) else unique(as.character(d_valid[[x_col]]))
  row_order <- if (is.factor(pvals_df[[y_col]])) levels(droplevels(pvals_df[[y_col]])) else unique(as.character(d_valid[[y_col]]))

  d_valid[[x_col]] <- factor(d_valid[[x_col]], levels = col_order)
  d_valid[[y_col]] <- factor(d_valid[[y_col]], levels = row_order)
  d_valid <- d_valid %>% dplyr::arrange(.data[[y_col]], .data[[x_col]], match(.data[[outcome_col]], outcome_levels))

  polys <- .build_split_circle_polygons(
    d = d_valid,
    x_col = x_col,
    y_col = y_col,
    id_cols = c(x_col, y_col),
    outcome_col = outcome_col,
    outcome_levels = outcome_levels,
    z_col = z_col,
    magnitude_col = magnitude_col,
    max_radius = max_radius
  )

  rings <- .build_significance_rings(
    d = d_valid,
    x_col = x_col,
    y_col = y_col,
    id_cols = c(x_col, y_col),
    magnitude_col = magnitude_col,
    difference_p_col = difference_p_col,
    p_threshold = difference_p_threshold,
    max_radius = max_radius
  )

  legend_pts <- d_valid %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(x_col, y_col)))) %>%
    dplyr::summarize(abs_mag = max(abs(.data[[magnitude_col]]), na.rm = TRUE), .groups = "drop")

  legend_pts$x_num <- match(as.character(legend_pts[[x_col]]), col_order)
  legend_pts$y_num <- match(as.character(legend_pts[[y_col]]), row_order)

  p <- ggplot2::ggplot()

  if (nrow(polys)) {
    p <- p + ggplot2::geom_polygon(
      data = polys,
      ggplot2::aes(x = x, y = y, group = poly_id, fill = z_value),
      linewidth = 0
    )
  }

  if (nrow(rings)) {
    p <- p + ggplot2::geom_polygon(
      data = rings,
      ggplot2::aes(x = x, y = y, group = poly_id),
      fill = NA,
      colour = "black",
      linewidth = 0.7
    )
  }

  p <- p +
    ggplot2::geom_point(
      data = legend_pts,
      ggplot2::aes(x = x_num, y = y_num, size = abs_mag),
      alpha = 0
    ) +
    ggplot2::scale_size_continuous(
      name = paste0( # "Max abs\n",
        magnitude_col),
      range = c(3, 8)
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(title = magnitude_lab,
                                   override.aes = list(alpha = 1, shape = 21, fill = "grey60", color = "gray60")
      )
    ) +
    ggplot2::scale_fill_gradient2(low = low,
                                  mid = mid,
                                  high = high,
                                  midpoint = 0,
                                  na.value = "grey85",
                                  name = fill_lab) +
    ggplot2::scale_x_continuous(breaks = seq_along(col_order),
                                labels = col_order,
                                position = "bottom",
                                expand = ggplot2::expansion(add = 0.3)) +
    ggplot2::scale_y_reverse(breaks = seq_along(row_order),
                             labels = row_order,
                             expand = ggplot2::expansion(add = 0.3)) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::labs(title = title,
                  x = x_lab,
                  y = y_lab,
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), plot.title = element_text(face = "plain"))

  # add custom half circle and significance ring legends
  main_legend <- cowplot::get_legend(
    p + ggplot2::theme(
      legend.position = "right",
      legend.box.spacing = ggplot2::unit(0, "pt"),
      legend.box = "vertical",
      legend.box.just = "left",
      legend.justification = "left",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      aspect.ratio = 1
    )
  )

  split_legend <- .make_half_circle_legend(
    outcome_levels = outcome_levels,
    title = split_lab,
    sig_label = paste0("FDR < ", difference_p_threshold)
  )

  legend_stack <- cowplot::plot_grid(
    main_legend,
    split_legend,
    ncol = 1,
    rel_heights = c(1, .32),
    align = "v",
    axis = "l"
  )

  cowplot::plot_grid(
    p + ggplot2::theme(legend.position = "none"),
    legend_stack,
    nrow = 1,
    rel_widths = c(1, .20),
    align = "h",
    axis = "tb"
  )

}

#' Plot spatial Tweedie contrasts as split circles
#'
#' Generates the split-circle visualization of band-vs-rest spatial contrasts
#' for each outcome group. Circle halves represent outcome-specific effects;
#' the black ring denotes a significant between-outcome difference.
#'
#' @param spatial_res Output from run_spatial_tweedie().
#' @param cfg Analysis configuration list.
#' @param outcome Outcome to plot; must occur in
#'   cfg$spatial_tweedie$outcome_levels.
#' @param save_plot Save the resulting PDF.
#'
#' @return Invisibly returns the prepared plotting data, plot, and output path.
plot_spatial_tweedie <- function(spatial_res, cfg,
                                 outcome = "Progression",
                                 save_plot = TRUE) {

  s <- cfg$spatial_tweedie

  if (is.null(s$outcome_levels[[outcome]]))
    stop("No outcome levels configured for `", outcome, "`.")
  if (!exists("plot_split_circle_grid", mode = "function"))
    stop("`plot_split_circle_grid()` must be defined before plotting.")

  outcome_levels <- s$outcome_levels[[outcome]]

  # ---- prepare final model contrasts for plotting ----------------------------
  d <- spatial_res$stats_long %>%
    dplyr::mutate(
      subtype_label = dplyr::case_when(
        subtype_label == "SP+/VGLUT1+ sc1" ~ "Cluster 1\n(SP+/VGLUT1+)",
        subtype_label == "TH.mid sc2" ~ "Cluster 2 (TH-mid)",
        subtype_label == "VACHT+ sc4" ~ "Cluster 4 (VACHT+)",
        subtype_label == "TH.high sc5" ~ "Cluster 5 (TH-high)",
        subtype_label == "SYP.high sc6" ~ "Cluster 6 (SYP-high)",
        subtype_label == "TH-/VACHT-/SP-/VGLUT1- sc7" ~
          "Cluster 7\n(Indeterminate)",
        subtype_label == "Vessel" ~ "Vessel (CD31+)",
        TRUE ~ as.character(subtype_label)
      ),
      band_label = gsub(" um", "", band_label),
      band_label = factor(
        band_label,
        levels = gsub(" um", "", s$band_levels)
      )
    )

  # ---- split-circle plot ------------------------------------------------------
  p <- plot_split_circle_grid(
    pvals_df = d,
    x_col = "band_label",
    y_col = "subtype_label",
    outcome_col = "outcome",
    outcome_levels = outcome_levels,

    # Tweedie log-link contrasts are natural-log ratios, NOT log2 fold changes.
    z_col = "band_vs_rest_log_ratio",
    magnitude_col = "outcome_log_ratio",

    difference_p_col = "padj_between_outcomes",
    difference_p_threshold = 0.05,

    title = paste(
      "Differential Spatial Enrichment in Relation to",
      gsub("_", " ", outcome)
    ),
    x_lab = "Binned distance from nearest tumor-stroma border (\u00b5m)",
    y_lab = NULL,
    fill_lab = "Log ratio\nEnrichment",
    magnitude_lab = "\u0394 log ratio\nEnrichment",
    split_lab = gsub("_", " ", outcome)
  )

  # ---- save ------------------------------------------------------------------
  out_pdf <- NULL
  if (save_plot) {
    out_dir <- file.path(cfg$path, cfg$paths$results$spatial, "figures")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    out_pdf <- file.path(
      out_dir,
      paste0(
        "fig7d_split_bubble_banded_nervegroup",
        ifelse(is.null(s$nerve_group_keep), "ALL", toupper(s$nerve_group_keep)),
        "_maxrad", s$max_radius_um, "_", outcome, ".pdf"
      )
    )

    ggplot2::ggsave(out_pdf, p, width = 8.25, height = 6)
    message("Saved split bubble plot to ", normalizePath(out_pdf, mustWork = FALSE))
  }

  invisible(list(
    data = d,
    plot = p,
    output_file = out_pdf
  ))
}


#' Plot raw STaN and vessel distances from the tumor-stroma border
#'
#' Uses the same images retained for the banded spatial analysis and plots the
#' distribution of individual STaN/vessel distances from the nearest tumor-
#' stroma border. Graphical settings are intentionally fixed.
#'
#' @param master_nerve_vessel Object-level nerve/vessel data.
#' @param spatial_res Output from run_spatial_tweedie().
#' @param cfg Analysis configuration list; only cfg$path is used here.
#' @param save_plots Save PDFs to the configured spatial figures directory.
#'
#' @return Invisibly returns the filtered plotting data and ggplot objects.
plot_spatial_distributions <- function(master_nerve_vessel, spatial_res, cfg,
                                       save_plots = TRUE, nerve_color_map) {

  # ---- plotting settings -----------------------------------------------
  xvar <- "dist_nearest_tumor_stroma_border_um"
  max_rad <- cfg$spatial_tweedie$max_radius_um
  outcomes <- c("BCG_failure", "Progression")

  subtype_levels <- rev(c(
    "Cluster 1 (SP+/VGLUT1+)", "Cluster 2 (TH-mid)",
    "Cluster 4 (VACHT+)", "Cluster 5 (TH-high)",
    "Cluster 6 (SYP-high)", "Cluster 7 (Indeterminate)",
    "Vessel (CD31+)"
  ))

  # Supply the manuscript palette explicitly; no replacement colors are inferred.
  pal <- nerve_color_map

  # ---- prepare object-level plotting data -----------------------------------
  d <- master_nerve_vessel %>%
    dplyr::left_join(
      spatial_res$spatial_clin %>%
        dplyr::select(image_id, BCG_failure, Progression, nerve_group),
      by = "image_id"
    ) %>%
    dplyr::filter(
      .data[[xvar]] <= max_rad,
      image_id %in% spatial_res$bands_keep$image_id,
      image_id %in% unique(spatial_res$nerve_border_bands$image_id),
      include_image_final_spatial_cohort %in% TRUE
    ) %>%
    dplyr::mutate(
      subtype_label = dplyr::case_when(
        subtype_label == "SP+/VGLUT1+ sc1" ~ "Cluster 1 (SP+/VGLUT1+)",
        subtype_label == "TH.mid sc2" ~ "Cluster 2 (TH-mid)",
        subtype_label == "VACHT+ sc4" ~ "Cluster 4 (VACHT+)",
        subtype_label == "TH.high sc5" ~ "Cluster 5 (TH-high)",
        subtype_label == "SYP.high sc6" ~ "Cluster 6 (SYP-high)",
        subtype_label == "TH-/VACHT-/SP-/VGLUT1- sc7" ~ "Cluster 7 (Indeterminate)",
        subtype_label == "Vessel" ~ "Vessel (CD31+)",
        TRUE ~ as.character(subtype_label)
      ),
      subtype_label = factor(subtype_label, levels = subtype_levels)
    )

  # ---- overall ridge plot ----------------------------------------------------
  p_overall <- ggplot2::ggplot(
    d, ggplot2::aes(
      x = .data[[xvar]], y = subtype_label,
      color = subtype_label, fill = subtype_label
    )
  ) +
    ggridges::stat_density_ridges(
      quantile_lines = TRUE, quantiles = 2, calc_ecdf = TRUE,
      jittered_points = TRUE,
      position = ggridges::position_points_jitter(width = 0.01, height = 0),
      point_shape = "|", point_size = 3, point_alpha = 1,
      alpha = 0.3, scale = 0.9, show.legend = FALSE, na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max_rad, 50), limits = c(0, max_rad),
      expand = c(0, 0)
    ) +
    cowplot::theme_cowplot() +
    ggplot2::labs(
      title = "Distance of STaNs and vessels to tumor-stroma border",
      x = "Distance to nearest tumor-stroma border (\u00b5m)",
      y = "Tissue type"
    )

  # ---- outcome-faceted ridge + density plots --------------------------------
  p_ridge <- p_density <- list()

  for (oc in outcomes) {
    dd <- d %>% dplyr::filter(!is.na(.data[[oc]]))

    p_ridge[[oc]] <- ggplot2::ggplot(
      dd, ggplot2::aes(
        x = .data[[xvar]], y = subtype_label,
        color = subtype_label, fill = subtype_label
      )
    ) +
      ggplot2::facet_wrap(stats::as.formula(paste0("~", oc)), ncol = 1) +
      ggridges::stat_density_ridges(
        quantile_lines = TRUE, quantiles = 2,
        jittered_points = TRUE,
        position = ggridges::position_points_jitter(width = 0.05, height = 0),
        point_shape = "|", point_size = 3, point_alpha = 1,
        alpha = 0.3, scale = 1, show.legend = FALSE, na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = pal) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::scale_x_continuous(
        breaks = seq(0, max_rad, 50), limits = c(0, max_rad)
      ) +
      cowplot::theme_cowplot() +
      ggplot2::labs(
        title = "Distance from STaNs and vessels to tumor-stroma border",
        x = "Distance to nearest tumor-stroma border (\u00b5m)",
        y = "Tissue type"
      )

    p_density[[oc]] <- ggpubr::ggdensity(
      dd, x = xvar, fill = "subtype_label", color = "subtype_label",
      alpha = 0, add = "median", rug = TRUE, palette = pal,
      facet.by = oc, ncol = 1, linewidth = 1, y = "density"
    ) +
      ggplot2::scale_x_continuous(
        breaks = seq(0, max_rad, 50), limits = c(0, max_rad)
      ) +
      cowplot::theme_cowplot() +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "plain")
      ) +
      ggplot2::labs(
        title = "Distance from STaNs and vessels to tumor-stroma border",
        x = "Distance to nearest tumor-stroma border (\u00b5m)",
        y = "Density"
      )
  }

  # ---- save ------------------------------------------------------------------
  if (save_plots) {
    out_dir <- file.path(cfg$path, cfg$paths$results$spatial, "figures")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    out_pdf <- file.path(
      out_dir, "fig7b_ridgeplot_dist_to_tumor.pdf"
    )
    ggplot2::ggsave(out_pdf, p_overall, width = 12, height = 7)
    message("Saved: ", normalizePath(out_pdf, mustWork = FALSE))

    for (oc in outcomes) {
      out_ridge <- file.path(
        out_dir, paste0("supp_spatial_ridge_", oc, ".pdf")
      )
      out_density <- file.path(
        out_dir, paste0("supp_spatial_density_", oc, ".pdf")
      )

      ggplot2::ggsave(out_ridge, p_ridge[[oc]], width = 8, height = 8)
      ggplot2::ggsave(out_density, p_density[[oc]], width = 8, height = 8)

      message("Saved distance ridge plot to ", normalizePath(out_ridge, mustWork = FALSE))
      message("Saved distance density plot to ", normalizePath(out_density, mustWork = FALSE))
    }
  }

  invisible(list(
    data = d,
    overall = p_overall,
    ridge_by_outcome = p_ridge,
    density_by_outcome = p_density
  ))
}
