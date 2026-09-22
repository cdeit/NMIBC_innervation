#' Plot histogram of STaN abundance or tissue area
#'
#' Generates a histogram showing the distribution of patient-level STaN+
#' area or quantified tissue area. Axis labels and plot titles are assigned
#' automatically for the standard tissue metrics.
#'
#' The resulting plot is saved as a PDF in the configured clinical figures directory
#' and returned as a ggplot object.
#'
#' @param df Data frame containing the variable to plot.
#' @param xvar Character string specifying the variable to plot. Supported
#'   variables include "nerve_obj_area_percent", "sum_stroma_area_um2",
#'   "sum_tumor_area_um2", and "sum_total_area_um2".
#' @param ylab Character string specifying the y-axis label.
#' @param title Optional plot title. For supported x variables, the title is
#'   assigned automatically.
#' @param binwidth Numeric histogram bin width.
#' @param fill_color Histogram fill color.
#' @param font_size Base font size used by cowplot.
#'
#' @return A ggplot object containing the histogram.
#'
#' @examples
#' plot_tissue_histogram(
#'   patient_nerve_quant,
#'   xvar = "nerve_obj_area_percent"
#' )
plot_tissue_histogram <- function(
    df,
    xvar = c("nerve_obj_area_percent", "sum_stroma_area_um2",
             "sum_tumor_area_um2", "sum_total_area_um2"),
    ylab = "Number of patients",
    title,
    binwidth = 0.025,
    font_size = 16, cfg) {

  # Set default variable
  if (is.null(xvar)) xvar <- "nerve_obj_area_percent"

  # Assign axis label and title based on the selected tissue metric
  if (xvar == "nerve_obj_area_percent") {
    xlab <- "STaN+ area (%)"
    title <- "Distribution of stromal tumor-associated nerve+ (STaN+) area"
    fill_color <- "orange"

  } else if (xvar == "sum_stroma_area_um2") {
    xlab <- "Stroma area (µm²)"
    title <- "Distribution of stromal area"
    fill_color <- "navy"

  } else if (xvar == "sum_tumor_area_um2") {
    xlab <- "Tumor area (µm²)"
    title <- "Distribution of tumor area"
    fill_color = "dodgerblue"

  } else if (xvar == "sum_total_area_um2") {
    xlab <- "Total tissue area (µm²)"
    title <- "Distribution of total tissue area"
    fill_color <- "darkgray"

  } else {
    xlab <- xvar
    title <- paste0("Distribution of ", xvar)
    fill_color <- "gray"
  }

  # Generate histogram
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[xvar]])) +
    ggplot2::geom_histogram(
      binwidth = binwidth,
      fill = fill_color,
      color = "white"
    ) +
    ggplot2::scale_x_continuous(expand = c(0.0075, 0.0075)) +
    ggplot2::scale_y_continuous(expand = c(0.025, 0.025)) +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      title = title
    ) +
    cowplot::theme_cowplot(font_size = font_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "plain")
    )

  # Create output directory and save histogram
  outdir <- file.path(cfg$path, cfg$paths$results$clinical, "figures")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  out_pdf <- file.path(outdir, paste0("histogram_", xvar, ".pdf"))
  ggplot2::ggsave(out_pdf, p, width = 7, height = 3.6)

  message("Saved histogram to: ", normalizePath(out_pdf))

  return(p)
}

#' Plot STaN abundance across clinical variables
#'
#' Creates the main clinical/outcome panel and supporting clinicopathologic panel.
#'
#' USAGE:
#'   df <- prep_patient_figure_data(patient_nerve_quant_07302025)
#'   res <- plot_stan_clinical_associations(df, cfg)
plot_stan_clinical_associations <- function(
    data, cfg,
    main_filename = "fig2_boxplots_stan_x_clinical.pdf",
    supp_filename = "supp_boxplots_fig_stan_x_clinical.pdf") {

  df <- data
  y <- "nerve_obj_area_percent"

  out_dir <- file.path(cfg$path, cfg$paths$results$clinical, "figures")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  main_pdf <- file.path(out_dir, main_filename)
  supp_pdf <- file.path(out_dir, supp_filename)

  paloc <- unname(c("gray90", "gray45"))
  pal2 <- unname(c("#0077B6", "#A4133C"))
  pal6 <- unname(c("#00A896", "#4361EE", "#A11D33", "#0B525B", "#9A6B2E", "#0077B6"))

  # Generic boxplot
  .box <- function(d, x, title, levels = NULL, labels = NULL) {
    d <- d %>% dplyr::filter(!is.na(.data[[x]]), !is.na(.data[[y]]))
    if (!is.null(levels)) d[[x]] <- factor(d[[x]], levels = levels)
    d <- d %>% dplyr::filter(!is.na(.data[[x]]))
    pal <- if (x %in% c("BCG_failure", "Progression")) paloc else if (dplyr::n_distinct(d[[x]]) <= 2) pal2 else pal6

    p <- ggplot2::ggplot(
      d, ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]])
    ) +
      ggplot2::geom_boxplot(outliers = FALSE, show.legend = FALSE) +
      ggplot2::geom_point(
        position = ggplot2::position_jitter(width = 0.2), show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = pal) +
      ggpubr::stat_compare_means(
        label = "p", show.legend = FALSE,
        method.args = list(exact = FALSE)
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(title = title, x = "", y = "STaN+ area (%)") +
      cowplot::theme_cowplot() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "plain"))

    if (!is.null(labels)) p <- p + ggplot2::scale_x_discrete(labels = labels)
    p
  }

  # Supporting clinicopathologic plots
  age <- ggplot2::ggplot(
    df %>% dplyr::filter(!is.na(Age)),
    ggplot2::aes(Age, nerve_obj_area_percent)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    # ggpmisc::stat_poly_line() +
    ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho") +
    ggplot2::labs(title = "Age at diagnosis", x = "Years", y = "STaN+ area (%)") +
    cowplot::theme_cowplot() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "plain"))


  sex <- .box(df %>% dplyr::filter(Sex != "Unknown"), "Sex", "Sex")
  smoking <- .box(df, "Smoking_group", "Smoking")
  cohort <- .box(df %>% dplyr::filter(!is.na(Cohort)), "Cohort", "RNA-seq cohort")
  retur <- .box(df, "reTURBT_prior_BCG", "re-TURBT prior to BCG induction", c("No", "Yes"))
  location <- .box(df, "Location_group", "Tumor location",
                   c("Wall", "Trigone/Base", "Dome", "Other"),
                   c("Wall" = "Bladder wall"))
  stage <- .box(df, "Tumor_stage", "Tumor stage", c("Ta", "T1"))
  size <- .box(df, "Size_group", "Tumor size", c("Small (<3cm)", "Large (>3cm)"))
  variant <- .box(df, "Variant_histology", "Variant histology")
  fgfr3 <- .box(df, "FGFR3_status", "FGFR3 mutation status", c("Wild-type", "Mutated"))
  lvi <- .box(df, "LVI", "Lymphovascular invasion")
  tils <- .box(df, "TILs_status", "Tumor-infiltrating lymphocytes")

  # Main clinical/outcome plots
  substage <- .box(df, "Substage", "T1 substage", c("T1m", "T1e"))
  cis <- .box(df, "Carcinoma_in_situ", "Carcinoma in situ", c("Absent", "Present"))
  focality <- .box(df, "Tumor_focality", "Tumor focality", c("Unifocal", "Multifocal"))
  bcg <- .box(df, "BCG_failure", "BCG failure")
  progression <- .box(df, "Progression", "Progression")

  # Criteria failure: pairwise Wilcoxon tests + BH correction
  crit_df <- df %>%
    dplyr::filter(
      Criteria_failure_clean %in% c("Relapsing", "Refractory", "MIBC"),
      !is.na(.data[[y]])
    ) %>%
    dplyr::mutate(
      Criteria_failure_clean = factor(
        Criteria_failure_clean,
        levels = c("Relapsing", "Refractory", "MIBC")
      )
    )

  criteria <- ggplot2::ggplot(
    crit_df,
    ggplot2::aes(
      x = Criteria_failure_clean,
      y = nerve_obj_area_percent,
      fill = Criteria_failure_clean
    )
  ) +
    ggplot2::geom_boxplot(outliers = FALSE, show.legend = FALSE) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.2), show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = pal6) +
    ggpubr::stat_pwc(method = "wilcox.test", p.adjust.method = "BH", na.rm = TRUE, exact = FALSE,
                     comparisons = list(c("Relapsing", "Refractory"),
                                        c("Relapsing", "MIBC"),
                                        c("Refractory", "MIBC"))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = "BCG failure criteria", x = "", y = "STaN+ area (%)") +
    cowplot::theme_cowplot() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "plain"))


  # Assemble
  main_patch <- (substage | cis | focality) / (bcg | progression | criteria)

  supp_patch <- (age | sex | smoking) /
    (cohort | retur | location) /
    (stage | size | variant) /
    (fgfr3 | lvi | tils)

  # Export
  ggplot2::ggsave(main_pdf, main_patch, width = 13, height = 8, units = "in")
  ggplot2::ggsave(supp_pdf, supp_patch, width = 13, height = 16, units = "in")

  message(
    "Saved:\n  ", normalizePath(main_pdf, mustWork = FALSE),
    "\n  ", normalizePath(supp_pdf, mustWork = FALSE)
  )

  invisible(list(
    main = main_patch, support = supp_patch,
    plots = list(
      age = age, sex = sex, smoking = smoking, cohort = cohort,
      reTURBT = retur, location = location, stage = stage, size = size,
      variant = variant, FGFR3 = fgfr3, LVI = lvi, TILs = tils,
      substage = substage, CIS = cis, focality = focality,
      BCG_failure = bcg, failure_criteria = criteria, progression = progression
    ),
    files = c(main = main_pdf, support = supp_pdf)
  ))
}

#' Plot STaN abundance across molecular subtypes and clinical outcomes
#'
#' Generates boxplots comparing STaN abundance across multiple molecular
#' subtype classifications and, within each subtype, across clinical outcomes.
#' Subtype columns, outcome recoding, filtering parameters, and color palettes
#' are supplied through `cfg`. Subtype levels with insufficient observations are
#' excluded automatically. Wilcoxon tests with BH correction are performed for
#' pairwise subtype comparisons and within-subtype outcome comparisons.
#'
#' The resulting subtype panels are combined into a single patchwork figure and
#' exported as a PDF. PDF height is fixed, while width scales automatically with
#' the total number of retained subtype groups.
#'
#' INPUTS:
#'   data     Data frame containing STaN abundance, subtype classifications,
#'            clinical outcomes, and filtering variables.
#'   cfg      Configuration list containing `subtype_plots` settings and
#'            `color_pals$subtypes` named palettes.
#'   filename Output PDF filename.
#'
#' OUTPUT:
#'   Invisibly returns a list containing:
#'     $plot  - combined patchwork figure
#'     $stats - list of statistical test results
#'     $data  - filtered/reformatted data used for plotting
#'     $pdf   - path to the exported PDF
#'   The PDF save location is also printed to the console.
#'
#' USAGE:
#'   res <- plot_nerves_by_subtype(patient_nerve_quant_07302025, cfg)
#'   res$plot
#'   res$stats

plot_nerves_by_subtype <- function(data, cfg, filename = "boxplots_nerves_x_subtypes_x_outcomes.pdf",
                                   outcome_colors = c("No" = "gray90", "Yes" = "gray45")) {
  # ----- settings/data setup -----
  s <- cfg$subtype_plots
  outdir <- file.path(cfg$path, cfg$paths$results$subtypes, "figures")
  base::dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # subtype_cols = named vector: plot name = source column
  subs <- names(s$subtype_cols)
  missing_cols <- base::setdiff(c(unname(s$subtype_cols), s$nerve_col, s$filter_col, s$followup_col, names(s$outcomes)), names(data))
  if (length(missing_cols)) base::stop("Missing columns: ", base::paste(missing_cols, collapse = ", "))

  df <- data %>%
    dplyr::filter(.data[[s$filter_col]] == s$filter_value, !is.na(.data[[s$followup_col]]))

  # duplicate source subtype columns under clean plotting names
  for (i in subs) df[[i]] <- df[[s$subtype_cols[[i]]]]

  # recode outcomes using named maps in cfg
  for (j in names(s$outcomes)) {
    map <- s$outcomes[[j]]
    df[[j]] <- dplyr::recode(base::as.character(df[[j]]), !!!map, .default = base::as.character(df[[j]]))
  }

  patches <- list(); patch_widths <- stats::setNames(numeric(length(subs)), subs); stats_list <- list()

  # ----- one 3-row patch per subtype system -----
  for (i in subs) {
    pal <- cfg$color_pals$subtypes[[i]]
    if (is.null(pal)) base::stop("No palette found at cfg$color_pals$subtypes$", i)

    # retain sufficiently represented levels; named palette controls level order when available
    counts <- table(df[[i]])
    levels_keep <- names(counts[counts > s$min_per_level])
    if (!is.null(names(pal))) levels_keep <- intersect(names(pal), levels_keep)

    df_ <- df %>%
      dplyr::filter(!is.na(.data[[i]]), .data[[i]] != "NA", .data[[i]] %in% levels_keep) %>%
      dplyr::mutate(!!i := factor(.data[[i]], levels = levels_keep))

    if (!nrow(df_)) next
    patch_widths[i] <- length(levels_keep)

    # ----- STaN abundance across subtype levels -----
    levs <- levels(df_[[i]])
    if (length(levs) > 1) {
      comps <- utils::combn(levs, 2, simplify = FALSE)
      stat_sub <- ggpubr::compare_means(
        data = df_, formula = stats::as.formula(paste(s$nerve_col, "~", i)),
        method = "wilcox.test", p.adjust.method = "BH", comparisons = comps
      ) %>%
        dplyr::mutate(
          comparison = paste0(group1, " vs ", group2),
          n_group1 = vapply(group1, \(x) sum(df_[[i]] == x), integer(1)),
          n_group2 = vapply(group2, \(x) sum(df_[[i]] == x), integer(1)),
          p.adj.method = "BH",
          p.adj.signif = dplyr::case_when(p.adj < .001 ~ "***", p.adj < .01 ~ "**", p.adj < .05 ~ "*", TRUE ~ "ns")
        ) %>%
        dplyr::relocate(comparison, .after = ".y.") %>%
        dplyr::relocate(n_group1, .after = group1) %>%
        dplyr::relocate(n_group2, .after = group2)

      stats_list[[paste0("nerves_x_", i)]] <- stat_sub
    }

    p1 <- ggplot2::ggplot(df_, ggplot2::aes(x = .data[[i]], y = .data[[s$nerve_col]], fill = .data[[i]])) +
      ggplot2::geom_boxplot(outliers = FALSE, show.legend = FALSE) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = .2), show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::labs(x = "Subtype", y = "STaN+ area (%)", title = gsub("_", " ", i)) +
      ggplot2::coord_cartesian(clip = "off") + cowplot::theme_cowplot()

    if (length(levs) > 2) {
      p1 <- p1 + ggpubr::stat_compare_means(show.legend = FALSE) +
        ggpubr::stat_pwc(method = "wilcox_test", p.adjust.method = "BH", hide.ns = TRUE,
                         label = "p.adj", position = ggplot2::position_jitter(width = .2), show.legend = NA)
    } else if (length(levs) == 2) p1 <- p1 + ggpubr::stat_compare_means(method = "wilcox.test", show.legend = FALSE)

    # ----- STaN abundance vs each clinical outcome, faceted by subtype -----
    outcome_plots <- list()
    for (j in names(s$outcomes)) {
      test_df <- df_ %>%
        dplyr::filter(!is.na(.data[[j]])) %>%
        dplyr::group_by(.data[[i]]) %>%
        dplyr::filter(dplyr::n_distinct(.data[[j]]) == 2,
                      min(table(.data[[j]])) >= s$min_per_outcome) %>%
        dplyr::ungroup()

      if (nrow(test_df)) {
        stat_df <- test_df %>%
          dplyr::rename(facet_var = dplyr::all_of(i), x_var = dplyr::all_of(j)) %>%
          dplyr::group_by(facet_var) %>%
          rstatix::wilcox_test(stats::as.formula(paste(s$nerve_col, "~ x_var"))) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(p.adj = p.adjust(p, "BH"), label = paste0("p.adj = ", signif(p.adj, 2)), outcome = j) %>%
          dplyr::relocate(outcome, .before = group1) %>%
          dplyr::rename(!!i := facet_var)

        stats_list[[paste0("nerves_x_", i, "_x_", j)]] <- stat_df
      }

      strip_cols <- unname(pal[levels_keep])
      p_ <- ggplot2::ggplot(df_, ggplot2::aes(x = .data[[j]], y = .data[[s$nerve_col]], fill = .data[[j]], group = .data[[j]])) +
        ggplot2::geom_boxplot(outliers = FALSE, show.legend = FALSE) +
        ggplot2::geom_point(position = ggplot2::position_jitter(width = .2), show.legend = FALSE) +
        ggplot2::scale_fill_manual(values = outcome_colors) +
        ggpubr::stat_compare_means(method = "wilcox.test", hide.ns = FALSE, label = "p",
                                   label.y.npc = .95) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.02, .15))) +
        ggplot2::labs(x = gsub("_", " ", j), y = "STaN+ area (%)") +
        cowplot::theme_cowplot() +
        ggh4x::facet_wrap2(stats::as.formula(paste("~", i)), nrow = 1,
                           strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = strip_cols, color = NA)))

      outcome_plots[[j]] <- p_
    }

    patch <- patchwork::wrap_plots(c(list(p1), outcome_plots), ncol = 1)
    patches[[i]] <- patch # patchwork::wrap_elements(patch)
    keep <- names(patches)

    quilt <- patchwork::wrap_plots(
      patches,
      nrow = 1,
      widths = patch_widths[keep]
    )
  }

  # ----- combine/export; width scales with total plotted subtype groups -----
  keep <- names(patches)
  quilt <- patchwork::wrap_plots(patches, nrow = 1, widths = patch_widths[keep])
  n_groups <- sum(patch_widths[keep])
  pdf_width <- max(30, 1.5 * n_groups)
  out_pdf <- file.path(outdir, filename)
  ggplot2::ggsave(out_pdf, quilt, width = pdf_width, height = 9.5, units = "in")
  base::message("Plot saved to: ", base::normalizePath(out_pdf))
  invisible(list(plot = quilt, stats = stats_list, data = df, pdf = out_pdf))

}



#' Plot STaN clinical and molecular subtype heatmaps
#'
#' Generates the STaN heatmaps used for Figure 2A, Figure 3A, and Figure S10C.
#' Patients are ordered by total STaN abundance. Annotation columns, cohort
#' filtering, cluster labels, and analytical thresholds are defined in
#' cfg$heatmaps.
#'
#' USAGE:
#'   res <- plot_stan_heatmaps(patient_nerve_quant_07302025, cfg, complete_cases = TRUE)
#'   res <- plot_stan_heatmaps(patient_nerve_quant_07302025, cfg, plots = "clinical", complete_cases = TRUE)
#'
plot_stan_heatmaps <- function(data, cfg, complete_cases = TRUE,
                               plots = c("clinical", "subtypes", "full")) {

  s <- cfg$heatmaps
  clinical_palettes <- list(general_outcomes = c("gray90", "gray45"), BCG_failure = c(No = "gray",
Yes = "gray10"), BCG = c(Responder = "gray", Failure = "gray10"
), Progression = c(`No Progression` = "gray", Progression = "gray10"
), Criteria_failure = c(Relapsing = "#009E73", Refractory = "#4D96FF",
MIBC = "#B5173A"), Cohort = c(A = "#1D3557", B = "#00A9A5", `Not sequenced` = "#E0E0E0"
), Sex = c(Female = "#A6BE54", Male = "#006400", Unknown = "#E0E0E0"
), Smoking = c(No = "#A6BE54", Yes = "#006400"), Re_TUR = c(No = "#A6BE54",
Yes = "#006400"), `reTURBT prior to BCG` = c(No = "#A6BE54",
Yes = "#006400"), Location = c(Other = "#8C7C6D", Wall = "#125a56",
`Trigone/Base` = "#F5B041", Dome = "#B23A48", Unknown = "#D8D8D8"
), Stage_grade = c(TaG3 = "#E7C6D7", T1G3 = "#8C2D5D"), `Stage/grade` = c(TaHG = "#E7C6D7",
T1HG = "#8C2D5D"), Substage = c(T1m = "#E7C6D7", T1e = "#8C2D5D"
), Size = c(Small = "#E17A5F", Large = "#A50026", Missing = "#E0E0E0"
), Tumor_size = c(Small = "#E17A5F", `Large > 3cm` = "#A50026",
Missing = "#E0E0E0"), Tumor_focality = c(Unifocal = "#E17A5F",
Multifocal = "#A50026"), CIS = c(NoCIS = "#E17A5F", Cis = "#A50026",
No = "#E17A5F", Yes = "#A50026"), Variant_histology = c(UCC = "#E17A5F",
`UCC + Variant` = "#A50026"), Variant = c(UCC = "#E17A5F", `UCC + Variant` = "#A50026"
), LVI = c(`No LVI` = "#E17A5F", LVI = "#A50026", No = "#E17A5F",
Yes = "#A50026"), TILs = c(No = "#E17A5F", Yes = "#A50026"),
    FGFR3_mut = c(No = "#E17A5F", Yes = "#A50026"))
  s$cluster_labels <- c(subcluster_1 = "Cluster 1 (VGLUT1+/SP+)", subcluster_2 = "Cluster 2 (TH-mid)",
subcluster_3 = "Cluster 3 (SYP+ non-neuronal)", subcluster_4 = "Cluster 4 (VAChT+)",
subcluster_5 = "Cluster 5 (TH-high)", subcluster_6 = "Cluster 6 (SYP-high)",
subcluster_7 = "Cluster 7 (Indeterminate/Other)")
  s$stan_group_colors <- c("STaN-low" = "#fff7bc", "STaN-medium" = "#fec44f", "STaN-high" = "#d95f0e")
  # Annotation palettes are intrinsic to this heatmap; see the definitions below.

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------

  # Rename selected input columns using names supplied in cfg.
  .select_rename <- function(df, cols) {
    out <- df[, unname(cols), drop = FALSE]
    names(out) <- names(cols)
    out
  }

  # Generate annotation colors, favoring palettes already defined in cfg.
  .ann_colors <- function(df, source_cols = NULL, palettes = list()) {
    out <- list()

    for (nm in names(df)) {
      src <- if (!is.null(source_cols) && nm %in% names(source_cols))
        unname(source_cols[nm]) else nm

      pal <- palettes[[nm]]
      if (is.null(pal)) pal <- palettes[[src]]

      if (!is.null(pal)) {
        out[[nm]] <- pal

      } else if (nm == "Age") {
        rng <- range(df$Age, na.rm = TRUE)
        out[[nm]] <- circlize::colorRamp2(c(rng[1], stats::median(rng), rng[2]),
                                          c("#F5F7F5", "#B5D1C8", "#14746F"))

      } else if (is.numeric(df[[nm]])) {
        x <- df[[nm]]
        q <- unique(as.numeric(stats::quantile(
          x, c(0, .25, .5, .75, 1), na.rm = TRUE
        )))

        if (length(q) == 1) q <- q + c(-1, 0, 1)

        cols <- grDevices::colorRampPalette(
          c("#1E5AA8", "#7FB3D5", "#ffffbf", "#D96B5F", "#A4133C")
        )(length(q))

        out[[nm]] <- circlize::colorRamp2(q, cols)

      } else {
        lv <- unique(stats::na.omit(as.character(df[[nm]])))
        out[[nm]] <- stats::setNames(
          grDevices::hcl.colors(length(lv), "Dark 3"), lv
        )
      }
      }

    out
  }

  # Concatenate ComplexHeatmap annotations vertically.
  .vcat <- function(...) {
    Reduce(function(x, y) ComplexHeatmap::`%v%`(x, y), list(...))
  }

  # ---------------------------------------------------------------------------
  # Prepare patient-level data
  # ---------------------------------------------------------------------------

  df <- data %>%
    dplyr::filter(
      .data[[s$stroma_area_col]] > s$min_stroma_area,
      .data[[s$filter_col]] == s$filter_value
    ) %>%

    # Standardize variables used only for plotting.
    dplyr::mutate(
      `STaN group` = dplyr::case_when(
        .data[[s$nerve_group_col]] == "Low" ~ "STaN-low",
        .data[[s$nerve_group_col]] == "High" ~ "STaN-high",
        is.na(.data[[s$nerve_group_col]]) ~ "STaN-medium",
        TRUE ~ as.character(.data[[s$nerve_group_col]])
      ),
      `STaN group` = factor(
        `STaN group`,
        levels = c("STaN-low", "STaN-medium", "STaN-high")
      ),

      CIS = dplyr::case_when(
        CIS == "CIS" | grepl("CIS", Stage_Grade) ~ "Yes",
        CIS == "NoCIS" | !grepl("CIS", Stage_Grade) ~ "No",
        TRUE ~ CIS
      ),

      Stage_Grade = dplyr::case_when(
        Stage_Grade %in% c("TaG3", "TaG3_CIS") ~ "TaHG",
        Stage_Grade %in% c("TaG2_HG", "TaG2_HG_CIS") ~ "TaHG",
        Stage_Grade %in% c("T1G2_HG", "T1G2_HG_CIS") ~ "T1HG",
        Stage_Grade %in% c("T1G3", "T1G3_CIS") ~ "T1HG",
        TRUE ~ Stage_Grade
      ),
      `Stage/grade` = Stage_Grade,

      Criteria_failure = dplyr::if_else(
        Criteria_failure == "Relapsing?", "Relapsing", Criteria_failure
      ),
      Cohort = dplyr::coalesce(as.character(Cohort), "Not sequenced"
      ),
      `STaN group` = dplyr::case_when(
        dplyr::ntile(.data[[s$nerve_col]], 3) == 1 ~ "STaN-low",
        dplyr::ntile(.data[[s$nerve_col]], 3) == 2 ~ "STaN-medium",
        dplyr::ntile(.data[[s$nerve_col]], 3) == 3 ~ "STaN-high"
      ),
      `STaN group` = factor(
        `STaN group`,
        levels = c("STaN-low", "STaN-medium", "STaN-high")
      ),
      Re_TUR = stringr::str_to_sentence(Re_TUR)
    ) %>%

    # Zero composition values in essentially nerve-negative tumors.
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with(s$composition_prefix),
        ~ dplyr::if_else(
          .data[[s$nerve_col]] < s$composition_zero_threshold,
          0, as.numeric(.x)
        )
      )
    )

  # Filter out patients with incomplete clinical annotations if requested
  if (complete_cases) {
    df <- df %>% dplyr::filter(!is.na(`RNA-seq`))
  } else {
    df <- df
  }

  # Order by STaN abundance for plotting
  df <- df %>%
    dplyr::arrange(.data[[s$nerve_col]]) %>%
    as.data.frame()

  sample_ids <- as.character(df[[s$patient_id_col]])

  # ---------------------------------------------------------------------------
  # Prepare annotation data
  # ---------------------------------------------------------------------------

  subtype_df <- .select_rename(df, s$subtype_cols)
  clin_df <- .select_rename(df, s$clin_cols)
  surv_df <- .select_rename(df, s$surv_cols)

  # STaN abundance.
  mat_total <- matrix(df[[s$nerve_col]], nrow = 1,
                      dimnames = list("STaN area", sample_ids))

  # STaN cluster composition
  prop_cols <- grep(
    paste0("^", s$composition_prefix),
    names(df),
    value = TRUE
  )

  prop_mat <- t(as.matrix(df[, prop_cols, drop = FALSE]))
  colnames(prop_mat) <- sample_ids

  # e.g. proportion_obj_in_subcluster_1 -> subcluster_1
  cluster_ids <- sub(
    paste0("^", s$composition_prefix),
    "",
    rownames(prop_mat)
  )

  # Get colors BEFORE replacing row names with display labels
  cluster_colors <- cfg$color_pals$nerve_subclusters[cluster_ids]

  # Convert internal IDs to publication labels
  cluster_labels <- s$cluster_labels[cluster_ids]
  cluster_labels[is.na(cluster_labels)] <- cluster_ids[is.na(cluster_labels)]

  rownames(prop_mat) <- unname(cluster_labels)
  names(cluster_colors) <- rownames(prop_mat)

  # Order by mean cluster abundance
  ord <- rev(order(rowMeans(prop_mat, na.rm = TRUE)))
  prop_mat <- prop_mat[ord, , drop = FALSE]
  cluster_colors <- cluster_colors[ord]

  # ---------------------------------------------------------------------------
  # Build annotation tracks
  # ---------------------------------------------------------------------------

  ha_surv <- ComplexHeatmap::HeatmapAnnotation(
    df = surv_df,
    col = .ann_colors(
      surv_df, s$surv_cols,
      clinical_palettes
    ),
    annotation_name_side = "left",
    na_col = "grey90",
    which = "column",
    show_annotation_name = TRUE
  )

  ha_clin <- ComplexHeatmap::HeatmapAnnotation(
    df = clin_df,
    col = .ann_colors(
      clin_df, s$clin_cols,
      clinical_palettes
    ),
    annotation_name_side = "left",
    na_col = "grey90",
    which = "column",
    show_annotation_name = TRUE
  )

  ha_subtype <- ComplexHeatmap::HeatmapAnnotation(
    df = subtype_df,
    col = .ann_colors(
      subtype_df, s$subtype_cols,
      cfg$color_pals$subtypes
    ),
    annotation_name_side = "left",
    na_col = "grey90",
    which = "column",
    show_annotation_name = TRUE
  )

  # Total STaN abundance bar.
  ha_total <- ComplexHeatmap::HeatmapAnnotation(
    `STaN+\narea (%)` = ComplexHeatmap::anno_barplot(
      t(mat_total),
      ylim = c(0, 1),
      baseline = 0,
      bar_width = 1,
      border = FALSE,
      gp = grid::gpar(fill = "orange", col = "white", lwd = 0),
      axis_param = list(
        side = "left",
        at = c(0, .25, .5, .75, 1),
        labels = c("0", "0.25", "0.5", "0.75", "1")
      )
    ),
    which = "column",
    annotation_name_side = "left",
    annotation_name_rot = 0,
    na_col = "gray90",
    height = grid::unit(20, "mm")
  )

  # STaN group + abundance + cluster composition.
  ha_bar <- ComplexHeatmap::HeatmapAnnotation(
    `STaN group` = df$`STaN group`,

    `STaN+\narea (%)` = ComplexHeatmap::anno_barplot(
      t(mat_total),
      ylim = c(0, 1),
      baseline = 0,
      bar_width = 1,
      border = FALSE,
      gp = grid::gpar(fill = "orange", col = "white"),
      axis_param = list(
        side = "left",
        at = c(0, .25, .5, .75, 1),
        labels = c("0", "0.25", "0.5", "0.75", "1")
      )
    ),

    `STaN\ncomposition` = ComplexHeatmap::anno_barplot(
      t(prop_mat),
      ylim = c(0, 1),
      baseline = 0,
      bar_width = 1,
      border = FALSE,
      gp = grid::gpar(fill = cluster_colors[rownames(prop_mat)], col = NA),
      axis_param = list(
        side = "left",
        at = c(0, .5, 1),
        labels = c("0", "0.5", "1")
      )
    ),

    which = "column",
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_rot = 0,
    show_legend = TRUE,
    na_col = "gray90",
    height = grid::unit(67, "mm"),
    gap = grid::unit(1.75, "mm"),
    col = list(`STaN group` = s$stan_group_colors)
  )

  leg_bar <- ComplexHeatmap::Legend(
    title = "STaN cluster",
    labels = names(cluster_colors),
    legend_gp = grid::gpar(fill = unname(cluster_colors)),
    direction = "vertical"
  )

  # ---------------------------------------------------------------------------
  # Assemble and export requested plots
  # ---------------------------------------------------------------------------

  plots <- match.arg(plots, c("clinical", "subtypes", "full"), several.ok = TRUE)

  specs <- list(
    clinical = list(
      ht = .vcat(ha_surv, ha_total, ha_clin),
      file = "fig2a_heatmap_stan_x_clin.pdf",
      width = 16, height = 6,
      ann_side = "right", heat_side = "right"
    ),
    subtypes = list(
      ht = .vcat(ha_surv, ha_total, ha_subtype),
      file = "fig3a_heatmap_stan_x_mol_subtypes.pdf",
      width = 15, height = 5,
      ann_side = "right", heat_side = "right"
    ),
    full = list(
      ht = .vcat(ha_surv, ha_bar, ha_clin, ha_subtype),
      file = "supp_fig10c_heatmap_clin_subtypes_stans_entire_imaging_cohort.pdf",
      width = 26, height = 9,
      ann_side = "left", heat_side = "right"
    )
  )

  out_files <- character()

  for (nm in plots) {
    z <- specs[[nm]]
    workflow <- if (nm == "subtypes") "subtypes" else "clinical"
    out_dir <- file.path(cfg$path, cfg$paths$results[[workflow]], "figures")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_pdf <- file.path(out_dir, z$file)

    grDevices::pdf(out_pdf, width = z$width, height = z$height, onefile = FALSE)

    ComplexHeatmap::draw(
      z$ht,
      annotation_legend_side = z$ann_side,
      heatmap_legend_side = z$heat_side,
      merge_legends = TRUE,
      heatmap_legend_list = if (nm == "full") list(leg_bar) else NULL
    )

    grDevices::dev.off()
    out_files[nm] <- out_pdf

    message("Saved: ", normalizePath(out_pdf, mustWork = FALSE))
  }

  invisible(list(
    data = df,
    heatmaps = lapply(specs[plots], `[[`, "ht"),
    files = out_files
  ))
}


#' Prepare patient-level variables used across manuscript figures
#'
#' Standardizes clinical labels used repeatedly by heatmaps, boxplots,
#' survival analyses, and other patient-level figures. Raw columns are retained.
#'
#' USAGE:
#'   df <- prep_patient_figure_data(patient_nerve_quant_07302025)
prep_patient_figure_data <- function(data) {

  data %>%
    dplyr::mutate(
      Criteria_failure_clean = dplyr::if_else(
        Criteria_failure == "Relapsing?", "Relapsing", as.character(Criteria_failure)
      ),

      # Tumor stage
      Tumor_stage = dplyr::case_when(
        grepl("^Ta", Stage_Grade) ~ "Ta",
        grepl("^T1", Stage_Grade) ~ "T1",
        TRUE ~ NA_character_
      ),
      Tumor_stage = factor(Tumor_stage, levels = c("Ta", "T1")),

      # Collapse stage/grade + CIS combinations for heatmap display
      Stage_grade = dplyr::case_when(
        Stage_Grade %in% c("TaG3", "TaG3_CIS") ~ "TaG3",
        Stage_Grade %in% c("TaG2_HG", "TaG2_HG_CIS") ~ "TaG2_HG",
        Stage_Grade %in% c("T1G2_HG", "T1G2_HG_CIS") ~ "T1G2_HG",
        Stage_Grade %in% c("T1G3", "T1G3_CIS") ~ "T1G3",
        TRUE ~ Stage_Grade
      ),

      # Reconcile explicit CIS field with CIS encoded in Stage_Grade
      Carcinoma_in_situ = if_else(grepl("CIS", Stage_Grade), "Present", "Absent"),
      Carcinoma_in_situ = factor(Carcinoma_in_situ, levels = c("Absent", "Present")
      ),

      # T1 substage
      Substage = factor(Substage, levels = c("T1m", "T1e")),

      # Smoking
      Smoking_group = dplyr::case_when(
        Smoking == "No" ~ "No/Stopped",
        Smoking == "Yes" ~ "Current",
        TRUE ~ as.character(Smoking)
      ),

      # TIL status
      TILs_status = dplyr::case_when(
        TILs == "Yes" ~ "Present",
        TILs == "No" ~ "Absent",
        TRUE ~ as.character(TILs)
      ),

      # FGFR3 status
      FGFR3_status = factor(
        dplyr::case_when(
          FGFR3_mut == "No" ~ "Wild-type",
          FGFR3_mut == "Yes" ~ "Mutated",
          TRUE ~ NA_character_
        ),
        levels = c("Wild-type", "Mutated")
      ),

      # Tumor focality
      Tumor_focality = factor(
        dplyr::na_if(as.character(Tumor_focality), "NA"),
        levels = c("Unifocal", "Multifocal")
      ),

      # Reconcile alternate tumor-size variables
      Size_group = factor(
        dplyr::case_when(
          Size == "Small" ~ "Small (<3cm)",
          (is.na(Size) | Size == "Missing") & Tumor_size == "Small" ~ "Small (<3cm)",
          Size == "Large" ~ "Large (>3cm)",
          (is.na(Size) | Size == "Missing") & Tumor_size == "Large >3cm" ~ "Large (>3cm)",
          TRUE ~ NA_character_
        ),
        levels = c("Small (<3cm)", "Large (>3cm)")
      ),

      # Consolidated tumor location
      Location_group = dplyr::case_when(
        grepl("wall", Location, ignore.case = TRUE) ~ "Wall",
        grepl("trigone|ureteral|neck|ureth", Location, ignore.case = TRUE) ~ "Trigone/Base",
        grepl("dome", Location, ignore.case = TRUE) ~ "Dome",
        is.na(Location) | Location == "Missing" ~ NA_character_,
        TRUE ~ "Other"
      ),

      # Standardize re-TURBT capitalization
      reTURBT_prior_BCG = dplyr::case_when(
        tolower(Re_TUR) == "yes" ~ "Yes",
        tolower(Re_TUR) == "no" ~ "No",
        TRUE ~ as.character(Re_TUR)
      ),

      # Preserve original Cohort but make a heatmap-friendly version
      RNAseq_cohort = dplyr::coalesce(as.character(Cohort), "Not sequenced"),

      # Variant histology display column
      Variant_histology = Variant,

      # BCG-failure subtype label
      Criteria_failure = dplyr::case_when(
        Criteria_failure == "Relapsing?" ~ "Relapsing",
        TRUE ~ as.character(Criteria_failure)
      )
    )
}
