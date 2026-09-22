# 07_compute_spatial_data.R
# Segmentation TIFFs -> single-channel masks -> spatial measurements -> validation.
# Usage from the repository root (automatically launches the Alpine container):
# Rscript scripts/07_compute_spatial_data.R GOALS PHENOTYPES [IMAGE_IDS] [RUN_NAME]
# Precomputed mode requires no container. Image processing requires an Alpine allocation.
source("config/config.R")
if (!is.logical(cfg$use_precomputed_spatial) ||
    length(cfg$use_precomputed_spatial) != 1L ||
    is.na(cfg$use_precomputed_spatial)) {
  stop("Set cfg$use_precomputed_spatial to TRUE or FALSE in config/config.R.")
}
if (cfg$use_precomputed_spatial) {
  message("Skipping spatial image processing: cfg$use_precomputed_spatial is TRUE. Workflow 08 uses its configured precomputed measurements.")
} else {
  args <- commandArgs(trailingOnly = TRUE)
  # Accept the previous explicit option for existing launch commands.
  if (length(args) && args[1] == "--process-images") args <- args[-1]
  if (!length(args) %in% 2:4) {
    stop("Usage: Rscript scripts/07_compute_spatial_data.R GOALS PHENOTYPES [IMAGE_IDS] [RUN_NAME]")
  }
  if (Sys.getenv("SPATIAL_ALPINE_WORKFLOW") != "1") {
    # Only dispatch here; image packages are loaded by the container child.
    status <- system2("bash", shQuote(c("scripts/run_spatial_from_raw_alpine.sh", args)))
    if (status != 0L) stop("Alpine spatial workflow failed; see launcher output.", call. = FALSE)
  } else {
    # Set by the launcher to avoid recursively starting another container.
    # The processing scripts still enforce their full runtime/package checks.
    source("R/spatial_image_workflow.R")
    run_spatial_image_processing(args, normalizePath("."))
  }
}
