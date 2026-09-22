# Orchestrate the verified image processing without changing measurement kernels.
# Run in the preserved environment; the shell launcher supplies container mounts.
run_spatial_image_processing <- function(args, root) {
  if (!length(args) %in% 2:4) {
    stop("Usage: 07_compute_spatial_data.R GOALS PHENOTYPES [IMAGE_IDS] [RUN_NAME]")
  }
  settings <- new.env(parent = baseenv())
  sys.source(file.path(root, "config/config.R"), envir = settings)
  paths <- settings$cfg$paths$spatial
  resolve <- function(path) {
    if (grepl("^(/|~)", path)) path.expand(path) else file.path(root, path)
  }
  goals <- normalizePath(args[1], mustWork = TRUE)
  phenotypes <- normalizePath(args[2], mustWork = TRUE)
  ids <- if (length(args) >= 3L) args[3] else resolve(paths$image_manifest)
  ids <- normalizePath(ids, mustWork = TRUE)
  name <- if (length(args) >= 4L) args[4] else "spatial_processing"
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", name)) stop("Invalid run name.")
  tissue <- resolve(paths$tissue_nerve_masks)
  vessel <- resolve(paths$vessel_masks)
  masks <- file.path(resolve(paths$generated_masks), name)
  run <- file.path(resolve(paths$measurement_runs), name)
  tables <- resolve(paths$measurements)
  if (!dir.exists(tissue) || !dir.exists(vessel)) stop("Missing segmentation input directory; see docs/spatial/README.md.")
  if (file.exists(masks) || file.exists(run)) stop("Choose a new run name; existing masks/runs are preserved.")
  schema <- new.env(parent = baseenv())
  sys.source(file.path(root, "R/output_schema.R"), envir = schema)
  filenames <- paste0(names(schema$output_schema), ".csv")
  if (!all(file.exists(file.path(goals, filenames)))) stop("Missing reference table.")
  if (any(file.exists(file.path(tables, filenames)))) stop("Preserve existing spatial measurement CSVs before publishing a replacement.")
  for (folder in c(dirname(masks), dirname(run), tables)) {
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  }
  run_step <- function(script, arguments) {
    status <- system2(file.path(R.home("bin"), "Rscript"),
      shQuote(c(file.path(root, script), arguments)))
    if (status != 0L) stop("Spatial image processing stopped at ", script, call. = FALSE)
  }
  run_step("R/prepare_spatial_masks.R", c(tissue, vessel, phenotypes, ids, masks))
  run_step("R/compute_spatial_measurements.R", c(masks, ids, run))
  run_step("R/verify_and_publish.R", c(run, goals, tables))
  # Preserve routing configuration alongside the existing scientific provenance.
  file.copy(file.path(root, "config/config.R"), file.path(run, "provenance", "repository_config.R"))
  message("Verified processed measurements saved to ", tables)
}
