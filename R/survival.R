#' ============================================================
# --- Survival analysis ---------------------------------------
#' ============================================================


#  --- HELPER FUNCTIONS --------------------------------------------------------

### --- Helpers for Kaplan-Meiers -----------------------
#' Prepare data for subtype-stratified Kaplan-Meier plots
#'
#' Internal helper used by both KM plotting functions. Applies cohort filters,
#' copies configured subtype columns, and generates median- or tertile-based
#' STaN abundance groups.
.prep_km_subtype_data <- function(data, cfg) {
  s <- cfg$subtype_km_plots
  if (is.null(s)) base::stop("Missing cfg$subtype_km_plots.")
  if (is.null(s$nerve_col)) base::stop("Missing cfg$subtype_km_plots$nerve_col.")
  if (is.null(s$subtype_cols)) base::stop("Missing cfg$subtype_km_plots$subtype_cols.")
  if (is.null(s$outcomes)) base::stop("Missing cfg$subtype_km_plots$outcomes.")

  strata <- if (is.null(s$strata)) "median" else s$strata
  s$min_per_level <- if (is.null(s$min_per_level)) 10 else s$min_per_level

  # confirm required columns exist
  outcome_cols <- unique(unlist(lapply(s$outcomes, \(x) c(x$time, x$event))))
  required <- unique(c(unname(s$subtype_cols), s$nerve_col, outcome_cols, names(s$filters)))
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols))
    base::stop("Missing columns in `data`: ", paste(missing_cols, collapse = ", "))

  # apply configured cohort filters
  df <- data
  for (x in names(s$filters))
    df <- df %>% dplyr::filter(.data[[x]] %in% s$filters[[x]])

  # copy source subtype columns to standardized plotting names
  for (j in names(s$subtype_cols))
    df[[j]] <- df[[s$subtype_cols[[j]]]]

  # use existing 2-level categorical strata, otherwise generate from numeric STaN abundance
  nerve_vals <- df[[s$nerve_col]]

  if ((is.character(nerve_vals) || is.factor(nerve_vals)) &&
     dplyr::n_distinct(nerve_vals, na.rm = TRUE) == 2) {

    df$nerve_content <- as.character(nerve_vals)
    lvls <- unique(stats::na.omit(df$nerve_content))

    # sensible ordering for common binary labels
    preferred <- c(
      "STaN-high", "STaN-low", "High STaN", "Low STaN",
      "High", "Low", "Above median", "Below median",
      "Positive", "Negative", "Yes", "No", "Pred STaN-high", "Pred STaN-low", "Predicted STaN-high", "Predicted STaN-low"
    )
    strata_levels <- c(intersect(preferred, lvls), setdiff(lvls, preferred))

  } else if (is.numeric(nerve_vals) && strata == "median") {
    med <- stats::median(nerve_vals, na.rm = TRUE)
    df$nerve_content <- ifelse(nerve_vals > med, "Above median", "Below median")
    strata_levels <- c("Above median", "Below median")

  } else if (is.numeric(nerve_vals) && strata == "tertile") {
    df$nerve_content <- dplyr::ntile(nerve_vals, 3)
    df$nerve_content <- dplyr::recode(
      as.character(df$nerve_content),
      `1` = "1st (low)", `2` = "2nd", `3` = "3rd (high)"
    )
    strata_levels <- c("1st (low)", "2nd", "3rd (high)")

  } else {
    base::stop(
      "`", s$nerve_col,
      "` must be numeric for median/tertile stratification or a 2-level character/factor."
    )
  }

  df$nerve_content <- factor(df$nerve_content, levels = strata_levels)

  list(data = df, settings = s, strata = strata, strata_levels = strata_levels)
}


### --- Helpers for CoxPH -------------------------------
.prep_cox_data <- function(data, cfg) {
  s <- cfg$nerve_subclusters

  # Estimate outcome associations within the study's adequate-BCG cohort,
  # retaining the established stromal-area eligibility threshold.
  df <- dplyr::filter(
    data,
    .data[[s$stroma_area_col]] > s$min_stroma_area,
    .data[[s$adequate_bcg_col]] == s$adequate_bcg_value
  )

  # Copy configured nerve predictors to their model/display names
  for (nm in names(s$cox_predictors))
    df[[nm]] <- df[[s$cox_predictors[[nm]]]]

  # Treat unknown sex as missing
  df$Sex[df$Sex == "Unknown"] <- NA

  # Set near-zero nerve values to zero and max-scale positive values
  for (v in names(s$cox_predictors)) {
    x <- as.numeric(df[[v]])
    x[x <= s$cox_zero_threshold] <- 0
    mx <- max(x[x > 0], na.rm = TRUE)
    df[[v]] <- if (is.finite(mx)) ifelse(x > 0, x / mx, 0) else 0
  }

  # Convert missing values to an explicit category where configured
  for (v in intersect(s$missing_as_level, names(df)))
    df[[v]] <- ifelse(
      is.na(df[[v]]), "Missing", as.character(df[[v]])
    )

  # Apply configured reference levels for categorical covariates
  for (v in names(s$factor_refs))
    if (v %in% names(df) &&
        s$factor_refs[[v]] %in% unique(df[[v]]))
      df[[v]] <- stats::relevel(
        factor(df[[v]]),
        s$factor_refs[[v]]
      )

  df
}

.tidy_cox <- function(fit, endpoint, model, keep = NULL) {
  z <- summary(fit)

  # Extract coefficients and calculate hazard ratios with 95% CIs
  out <- tibble::tibble(
    term = rownames(z$coefficients),
    beta = z$coefficients[, "coef"],
    se = z$coefficients[, "se(coef)"],
    p = z$coefficients[, "Pr(>|z|)"],
    endpoint,
    model
  ) %>%
    dplyr::mutate(
      HR = exp(beta),
      lo = exp(beta - 1.96 * se),
      hi = exp(beta + 1.96 * se)
    ) %>%
    dplyr::select(endpoint, model, term, HR, lo, hi, p)

  # Optionally retain only the predictor of interest
  if (!is.null(keep))
    dplyr::filter(out, term == keep)
  else
    out
}

.fit_cox_set <- function(df, vars, time, event, endpoint, covars = NULL) {

  # Fit one Cox model per configured nerve predictor
  purrr::map_dfr(vars, function(v) {

    # Restrict each model to complete cases
    use <- c(time, event, v, covars)
    dat <- tidyr::drop_na(
      dplyr::select(df, dplyr::all_of(use))
    )

    # Skip models with insufficient observations or events
    if (nrow(dat) < 5 || sum(dat[[event]]) < 2)
      return(NULL)

    # Construct model formula from predictor and optional covariates
    rhs <- paste(sprintf("`%s`", c(v, covars)), collapse = " + ")
    form <- stats::as.formula(
      sprintf(
        "survival::Surv(`%s`, `%s`) ~ %s",
        time, event, rhs
      )
    )

    fit <- survival::coxph(form, data = dat)

    # Extract predictor-level model statistics
    .tidy_cox(
      fit,
      endpoint,
      if (length(covars)) "multivariable" else "univariate",
      if (length(covars)) paste0("`", v, "`") else NULL
    )
  }) %>%

    # Correct predictor P values within each model set
    dplyr::mutate(
      FDR = p.adjust(p, method = "BH"),
      term = gsub("`", "", term, fixed = TRUE)
    ) %>%
    dplyr::arrange(FDR)
}

.forest_plot <- function(res, title, subtitle, pal) {

  # Preserve result ordering from the model table
  res <- dplyr::mutate(
    res,
    term = factor(term, levels = rev(unique(term)))
  )

  # Plot hazard ratios and 95% confidence intervals
  ggplot2::ggplot(
    res,
    ggplot2::aes(HR, term, fill = term)
  ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = 2
    ) +
    ggplot2::geom_errorbarh(
      orientation = "y",
      ggplot2::aes(xmin = lo, xmax = hi),
      height = 0.2
    ) +
    ggplot2::geom_point(
      size = 6,
      shape = 21
    ) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Hazard ratio (log10 scale)",
      y = NULL
    ) +
    cowplot::theme_cowplot() +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      plot.margin = ggplot2::margin(10, 120, 10, 10),
      legend.position = "none"
    ) +

    # Add HR and FDR values beside each estimate
    ggplot2::geom_text(
      ggplot2::aes(
        x = Inf,
        label = sprintf(
          "HR = %.2f   (%s)", HR,
          ifelse(FDR < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", FDR)))
      ),
      hjust = -0.1,
      size = 3.2
    )
}

.run_cox_endpoint <- function(df, time, event, label, s, pal) {
  vars <- names(s$cox_predictors)

  # Fit univariate and covariate-adjusted models
  uni <- .fit_cox_set(
    df, vars, time, event, label)

  multi <- .fit_cox_set(df, vars, time, event, label, s$cox_adjust_covariates)

  # Return model tables and corresponding forest plots
  list(
    univariate_table = uni,
    multivariable_table = multi,

    univariate_forest = .forest_plot(uni, label, "Univariate Cox PH", pal
    ),

    multivariable_forest = .forest_plot(multi, label,
                                        paste("Multivariable Cox PH adjusted for",
                                              paste(s$cox_adjust_covariates, collapse = ", ")),
                                        pal
    ))
}



## --- KAPLAN-MEIER PLOTS -------------------------------------------------------

#' Plot basic STaN-stratified Kaplan-Meier curves
#'
#' Generates Kaplan-Meier curves for the configured survival outcomes using a
#' single grouping variable. Group levels are inferred automatically using a
#' preferred ordering for common binary labels. Log-rank results and median
#' survival estimates are printed to the console.
#'
#' @param data Data frame containing survival and grouping variables.
#' @param group_col Character string specifying the grouping variable.
#' @param cfg Analysis configuration list containing survival outcome settings.
#' @param palette Colors corresponding to the inferred group levels.
#' @param filename Output PDF filename.
#' @param workflow Result-root key for the generic KM plot (default clinical).
#' @param break_x X-axis break interval in months.
#' @param xlim X-axis limits.
#'
#' @return Invisibly returns a list of ggsurvplot objects, one per outcome.
plot_km_basic <- function(data, group_col, cfg,
                          palette = cfg$color_pals$nerve_strata,
                          filename = "KM_basic.pdf",
                          break_x = 24, xlim = c(0, 72), workflow = "clinical") {

  # Infer sensible group ordering
  lvls <- unique(as.character(data[[group_col]]))
  preferred <- c(
    "STaN-high", "STaN-low", "High STaN", "Low STaN",
    "High", "Low", "Above median", "Below median",
    "Positive", "Negative", "Yes", "No",
    "Pred STaN-high", "Pred STaN-low",
    "Predicted STaN-high", "Predicted STaN-low"
  )
  strata_levels <- c(intersect(preferred, lvls), setdiff(lvls, preferred))

  plots <- lapply(names(cfg$survival_outcomes), function(outcome) {

    s <- cfg$survival_outcomes[[outcome]]

    # Prepare complete survival data
    df <- data.frame(
      Time = data[[s$time]],
      Event = data[[s$event]],
      Group = factor(data[[group_col]], levels = strata_levels)
    ) |>
      dplyr::filter(!is.na(Time), !is.na(Event), !is.na(Group))

    # Fit KM curves and log-rank test
    km_fit <- survival::survfit(
      survival::Surv(Time, Event) ~ Group,
      data = df
    )

    surv_diff <- survival::survdiff(
      survival::Surv(Time, Event) ~ Group,
      data = df
    )

    p_value <- 1 - stats::pchisq(
      surv_diff$chisq,
      length(surv_diff$n) - 1
    )

    # Print survival statistics
    cat("\n", outcome, "\n", strrep("-", 40), "\n", sep = "")
    print(surv_diff)

    median_surv <- summary(km_fit)$table[
      , c("records", "events", "median", "0.95LCL", "0.95UCL"),
      drop = FALSE
    ]
    print(median_surv)
    cat("Log-rank p =", signif(p_value, 4), "\n")

    # Generate KM plot
    survminer::ggsurvplot(
      km_fit,
      data = df,
      title = s$title,
      pval = TRUE,
      conf.int = TRUE,
      palette = palette,
      break.x.by = break_x,
      xlab = "Time (months)",
      xlim = xlim,
      risk.table = TRUE,
      risk.table.y.text = FALSE,
      risk.table.height = 0.25,
      legend.labs = strata_levels,
      legend.title = "",
      legend = "top"
    )
  })

  names(plots) <- names(cfg$survival_outcomes)

  # Arrange KM plots and risk tables side-by-side
  combined <- survminer::arrange_ggsurvplots(
    plots,
    ncol = length(plots),
    nrow = 1,
    print = FALSE
  )

  # Save PDF
  outdir <- file.path(cfg$path, cfg$paths$results[[workflow]], "figures")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  out_pdf <- file.path(outdir, filename)

  ggplot2::ggsave(
    out_pdf,
    combined,
    width = 7 * length(plots),
    height = 6
  )

  message("Saved KM plot to: ", normalizePath(out_pdf))

  print(combined)
  invisible(plots)
}


#' Plot subtype + STaN-stratified Kaplan-Meier curves
#'
#' Generates one PDF page per molecular subtype system, with RFS and PFS shown
#' side-by-side. Curves are colored by subtype and linetyped by STaN abundance.
#'
#' USAGE:
#'   res <- plot_km_subtype_overlay(patient_nerve_quant_07302025, cfg)
plot_km_subtype_overlay <- function(data, cfg,
                                    filename = "KM_subtype_x_STaN_overlay.pdf", xlim = c(0, 72)) {

  prep <- .prep_km_subtype_data(data, cfg)
  df <- prep$data; s <- prep$settings

  outdir <- file.path(cfg$path, cfg$paths$results$subtypes, "figures")
  base::dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  out_pdf <- file.path(outdir, filename)

  plots <- list(); fits <- list()

  # open multipage PDF
  # if file exists, delete it first
  if (file.exists(out_pdf)) base::file.remove(out_pdf)
  grDevices::pdf(out_pdf, # width = 14, height = 6.5
                 width = 20, height = 12,
                 onefile = TRUE)

  for (j in names(s$subtype_cols)) {

    subtype_pal <- cfg$color_pals$subtypes[[j]]
    if (is.null(subtype_pal))
      base::stop("Missing cfg$color_pals$subtypes$", j)

    # retain sufficiently represented subtype levels
    counts <- table(df[[j]])
    levels_keep <- names(counts[counts > s$min_per_level])

    if (!is.null(names(subtype_pal)))
      levels_keep <- intersect(names(subtype_pal), levels_keep)

    if (!length(levels_keep)) next

    df_ <- df %>%
      dplyr::filter(
        !is.na(.data[[j]]),
        .data[[j]] != "NA",
        .data[[j]] %in% levels_keep
      ) %>%
      dplyr::mutate(
        Subtype = factor(.data[[j]], levels = levels_keep)
      )

    subtype_cols <- unname(subtype_pal[levels_keep])
    outcome_plots <- list()

    for (i in names(s$outcomes)) {

      o <- s$outcomes[[i]]

      df_surv <- df_ %>%
        dplyr::filter(
          !is.na(.data[[o$time]]),
          !is.na(.data[[o$event]]),
          !is.na(nerve_content)
        ) %>%
        dplyr::mutate(
          .time = .data[[o$time]],
          .event = if (!is.null(o$event_map))
            unname(o$event_map[as.character(.data[[o$event]])])
          else as.numeric(.data[[o$event]])
        ) %>%
        dplyr::filter(!is.na(.event))

      fit <- survival::survfit(
        survival::Surv(.time, .event) ~ Subtype + nerve_content, data = df_surv)

      p <- survminer::ggsurvplot(
        fit,
        data = df_surv,
        color = "Subtype",
        palette = subtype_cols,
        linetype = "nerve_content",
        pval = TRUE,
        xlab = "Time (months)",
        ylab = i,
        title = gsub("_", " ", j),
        xlim = xlim,
        break.time.by = 24,
        risk.table = TRUE,
        risk.table.y.text = TRUE,
        risk.table.y.text.col = FALSE,
        risk.table.height = 0.3
      )

      fits[[j]][[i]] <- fit
      plots[[j]][[i]] <- p
      outcome_plots[[i]] <- p
    }

    # print one page: RFS + PFS side-by-side
    page <- survminer::arrange_ggsurvplots(
      outcome_plots,
      nrow = 1,
      ncol = length(outcome_plots),
      print = FALSE
    )

    print(page)
  }

  grDevices::dev.off()

  base::message(
    "Plot saved to: ", base::normalizePath(out_pdf)
  )

  invisible(list(
    plots = plots,
    fits = fits,
    data = df,
    pdf = out_pdf
  ))
}


#' Plot STaN-stratified KM curves faceted by molecular subtype
#'
#' Generates one PDF page per molecular subtype system. Each page contains the
#' configured survival outcomes side-by-side, with individual subtype groups
#' shown as facets and subtype-colored ggh4x strip labels. Curves are colored
#' by STaN abundance stratum.
#'
#' USAGE:
#'   res <- plot_km_subtype_faceted(patient_nerve_quant_07302025, cfg)
plot_km_subtype_faceted <- function(data, cfg,
                                    filename = "KM_STaN_faceted_by_subtype.pdf", xlim = c(0, 72)) {

  prep <- .prep_km_subtype_data(data, cfg)
  df <- prep$data; s <- prep$settings

  # STaN palette in explicit factor order
  strata_pal <- cfg$color_pals$nerve_strata[[prep$strata]]
  if (is.null(strata_pal))
    base::stop("Missing cfg$color_pals$nerve_strata$", prep$strata)

  if (!is.null(names(strata_pal)))
    strata_pal <- unname(strata_pal[prep$strata_levels])

  outdir <- file.path(cfg$path, cfg$paths$results$subtypes, "figures")
  base::dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  out_pdf <- file.path(outdir, filename)

  plots <- list(); fits <- list()

  # determine retained subtype groups for each system before opening PDF
  retained_levels <- lapply(names(s$subtype_cols), function(j) {
    subtype_pal <- cfg$color_pals$subtypes[[j]]
    if (is.null(subtype_pal))
      base::stop("Missing cfg$color_pals$subtypes$", j)

    counts <- table(df[[j]])
    levels_keep <- names(counts[counts > s$min_per_level])

    # palette order determines facet order
    if (!is.null(names(subtype_pal)))
      levels_keep <- intersect(names(subtype_pal), levels_keep)

    levels_keep
  })
  names(retained_levels) <- names(s$subtype_cols)

  n_groups <- lengths(retained_levels)
  max_groups <- max(n_groups, na.rm = TRUE)

  if (max_groups < 1)
    base::stop("No molecular subtype groups passed the minimum sample-size filter.")

  # fixed physical width per KM facet; PDF width scales to largest subtype system
  facet_width <- 2
  n_outcomes <- length(s$outcomes)
  pdf_width <- 2 + (n_outcomes * max_groups * facet_width)
  pdf_height <- 2.75

  if (grDevices::dev.cur() > 1) grDevices::graphics.off()
  if (file.exists(out_pdf)) base::file.remove(out_pdf)

  grDevices::pdf(
    out_pdf,
    width = pdf_width,
    height = pdf_height,
    onefile = TRUE
  )
  on.exit(if (grDevices::dev.cur() > 1) grDevices::dev.off(), add = TRUE)

  first_page <- TRUE

  # one molecular subtype system per PDF page
  for (j in names(s$subtype_cols)) {

    subtype_pal <- cfg$color_pals$subtypes[[j]]
    levels_keep <- retained_levels[[j]]

    if (!length(levels_keep)) next

    df_ <- df %>%
      dplyr::filter(
        !is.na(.data[[j]]),
        .data[[j]] != "NA",
        .data[[j]] %in% levels_keep
      ) %>%
      dplyr::mutate(
        Subtype = factor(.data[[j]], levels = levels_keep),
        nerve_content = factor(nerve_content, levels = prep$strata_levels)
      )

    strip_cols <- unname(subtype_pal[levels_keep])
    outcome_plots <- list()

    # RFS/PFS
    for (i in names(s$outcomes)) {

      o <- s$outcomes[[i]]

      df_surv <- df_ %>%
        dplyr::filter(
          !is.na(.data[[o$time]]),
          !is.na(.data[[o$event]]),
          !is.na(nerve_content)
        ) %>%
        dplyr::mutate(
          .time = .data[[o$time]],
          .event = if (!is.null(o$event_map))
            unname(o$event_map[as.character(.data[[o$event]])])
          else as.numeric(.data[[o$event]])
        ) %>%
        dplyr::filter(!is.na(.event))

      fit <- survival::survfit(
        survival::Surv(.time, .event) ~ Subtype + nerve_content,
        data = df_surv
      )

      # faceted survminer plot when >1 retained subtype
      if (length(levels_keep) > 1) {

        p <- survminer::ggsurvplot_facet(
          fit,
          data = df_surv,
          facet.by = "Subtype",
          color = "nerve_content",
          palette = strata_pal,
          pval = TRUE,
          risk.table = TRUE,
          risk.table.y.text = FALSE,
          xlab = "Time (months)",
          ylab = i,
          title = gsub("_", " ", j),
          xlim = xlim,
          break.time.by = s$break_time,
          ggtheme = ggplot2::theme_classic()
        ) +
          ggh4x::facet_wrap2(
            ~Subtype,
            nrow = 1,
            strip = ggh4x::strip_themed(
              background_x = ggh4x::elem_list_rect(
                fill = strip_cols,
                color = strip_cols
              )
            )
          )

      } else {

        # standard survminer plot for one retained subtype
        g <- survminer::ggsurvplot(
          fit,
          data = df_surv,
          color = "nerve_content",
          palette = strata_pal,
          pval = TRUE,
          xlab = "Time (months)",
          ylab = i,
          title = gsub("_", " ", j),
          xlim = xlim,
          break.time.by = s$break_time,
          surv.median.line = "hv",
          ggtheme = ggplot2::theme_classic()
        )

        p <- g$plot +
          ggh4x::facet_wrap2(
            ~Subtype,
            nrow = 1,
            strip = ggh4x::strip_themed(
              background_x = ggh4x::elem_list_rect(
                fill = strip_cols,
                color = strip_cols
              )
            )
          )
      }

      fits[[j]][[i]] <- fit
      plots[[j]][[i]] <- p

      # pad narrower subtype systems so individual facet width stays constant
      n_current <- length(levels_keep)
      pad <- max_groups - n_current

      if (pad > 0) {
        p <- cowplot::plot_grid(
          NULL, p, NULL,
          nrow = 1,
          rel_widths = c(pad / 2, n_current, pad / 2)
        )
      }

      outcome_plots[[i]] <- p
    }

    # configured outcomes side-by-side on one page
    page <- cowplot::plot_grid(
      plotlist = outcome_plots,
      nrow = 1,
      align = "hv"
    )

    if (!first_page) grid::grid.newpage()

    grid::grid.draw(
      ggplot2::ggplotGrob(page)
    )

    first_page <- FALSE
  }

  grDevices::dev.off()

  base::message(
    "Plot saved to: ", base::normalizePath(out_pdf),
    "\nPDF dimensions: ", round(pdf_width, 1), " x ", pdf_height, " in per page",
    "\nMaximum retained subtype groups: ", max_groups
  )

  invisible(list(
    plots = plots,
    fits = fits,
    data = df,
    pdf = out_pdf
  ))
}


## --- COX PROPORTIONAL HAZARDS MODELING ----------------------------------------

#' Run nerve-subcluster Cox proportional hazards models
#'
#' Fits univariate and covariate-adjusted Cox models for each configured nerve
#' predictor across all configured survival endpoints. P values are BH-adjusted
#' within each endpoint/model set. Forest plots and model tables are exported.
#'
#' @param data Patient-level data containing survival outcomes, nerve predictors,
#'   cohort-filtering variables, and adjustment covariates.
#' @param cfg Analysis configuration; see `config.R`.
#' @param filename Multipage forest-plot PDF filename.
#' @param table_filename Excel filename containing model results.
#'
#' @return Invisibly returns a named list of model results for each endpoint.
#'
#' @examples
#' res_cox <- run_cluster_cox(patient_nerve_quant_07302025, cfg)
run_cluster_cox <- function(data, cfg,
                            filename = "subcluster_cox_forest_plots.pdf",
                            table_filename = "subcluster_cox_results.xlsx") {

  s <- cfg$nerve_subclusters
  df <- .prep_cox_data(data, cfg)

  # Fit univariate and multivariable models for each survival endpoint
  results <- lapply(names(s$cox_endpoints), function(endpoint) {
    ep <- s$cox_endpoints[[endpoint]]
    .run_cox_endpoint(
      df, ep$time, ep$event, endpoint, s,
      cfg$color_pals$nerve_subclusters
    )
  })
  names(results) <- names(s$cox_endpoints)

  # Collect and save forest plots
  plots <- unlist(
    lapply(results, function(x)
      list(x$univariate_forest, x$multivariable_forest)),
    recursive = FALSE
  )

  out_pdf <- file.path(cfg$path, cfg$paths$results$nerve_clusters, "figures", filename)
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

  pages <- gridExtra::marrangeGrob(plots, nrow = 1, ncol = 1)
  ggplot2::ggsave(out_pdf, pages, width = 9, height = 4, units = "in")

  # Collect model tables and save one sheet per endpoint/model type
  sheets <- list()
  for (nm in names(results)) {
    sheets[[paste(nm, "univariate")]] <- results[[nm]]$univariate_table
    sheets[[paste(nm, "multivariable")]] <- results[[nm]]$multivariable_table
  }

  out_xlsx <- file.path(cfg$path, cfg$paths$results$nerve_clusters, "tables", table_filename)
  dir.create(dirname(out_xlsx), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(sheets, out_xlsx)

  message("Saved: ", normalizePath(out_pdf, mustWork = FALSE))
  message("Saved: ", normalizePath(out_xlsx, mustWork = FALSE))

  invisible(results)
}
