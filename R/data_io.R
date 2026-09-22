# Shared input checks for the manuscript workflows.

#' Require the files supplied to a workflow
require_input_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Required input files are missing:\n  ", paste(missing, collapse = "\n  "),
      "\nSee README.md and docs/REPRODUCIBILITY.md for input provenance.",
      call. = FALSE
    )
  }
  invisible(paths)
}

#' Read a named object from an archived R workspace
#'
#' Several distributed files have an .rds suffix but were written with save().
#' Load them into an isolated environment; do not change the supplied datasets.
read_workspace_object <- function(path, object_name) {
  require_input_files(path)
  data_env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = data_env)
  if (!object_name %in% object_names) {
    stop("Expected object '", object_name, "' in ", path, call. = FALSE)
  }
  get(object_name, envir = data_env, inherits = FALSE)
}
