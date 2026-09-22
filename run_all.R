# Run selected manuscript workflows in separate R sessions.
# Run from the repository root:
# Rscript --vanilla run_all.R 02 03 04 06
# Steps run in the supplied order; failure stops execution.
# Select upstream regeneration (01/07) explicitly. Inputs: data/README.md.
# Invoke a script directly for its optional flags (see README.md).
steps <- commandArgs(trailingOnly = TRUE)
if (!length(steps) || any(!steps %in% sprintf("%02d", 1:8))) {
  stop("Select workflows explicitly, e.g. Rscript run_all.R 02 04 06. See README.md.", call. = FALSE)
}
for (step in steps) {
  script <- list.files("scripts", pattern = paste0("^", step, "_.*\\.R$"), full.names = TRUE)
  stopifnot(length(script) == 1L)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)))
  if (status != 0L) stop("Workflow failed: ", script, call. = FALSE)
}
