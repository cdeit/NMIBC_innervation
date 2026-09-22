## io.R -- mask I/O, clinical/manual_eval loading, image<->patient metadata.

## All masks in this pipeline are written with foreground = -1 (binary) or the
## negated subtype code (coded). Verified against real files this session --
## the assertions below make a future upstream encoding change fail loudly
## instead of silently corrupting every downstream measurement.

enr_read_binary <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  m <- tiff::readTIFF(path, as.is = TRUE)
  if (length(dim(m)) == 3) m <- m[, , 1]
  storage.mode(m) <- "integer"
  vals <- sort(unique(as.vector(m)))
  if (!all(vals %in% c(0L, -1L))) {
    stop("Unexpected encoding in binary mask ", path, ": values = ",
         paste(vals, collapse = ","))
  }
  (-m) > 0
}

enr_read_coded <- function(path, valid_codes = c(0L, 1L, 2L, 4L, 5L, 6L, 7L)) {
  if (!file.exists(path)) stop("File not found: ", path)
  m <- tiff::readTIFF(path, as.is = TRUE)
  if (length(dim(m)) == 3) m <- m[, , 1]
  storage.mode(m) <- "integer"
  coded <- -m
  vals <- sort(unique(as.vector(coded)))
  if (!all(vals %in% valid_codes)) {
    stop("Unexpected encoding in coded mask ", path, ": values = ",
         paste(vals, collapse = ","))
  }
  coded
}

# Only five single-channel inputs are required.
layer_paths <- function(image_id, cfg) {
  d <- file.path(cfg$in_dir, image_id)
  list(tumor=file.path(d,cfg$file_tumor), stroma=file.path(d,cfg$file_stroma),
       glass=file.path(d,cfg$file_glass), vessel=file.path(d,cfg$file_vessel),
       phenotyped=file.path(d,cfg$file_nerve_phenotyped))
}
