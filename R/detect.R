.sc_mtx_stem <- function(path) {
  nm <- tolower(basename(path))
  nm <- sub("\\.gz$", "", nm)
  sub("\\.mtx$", "", nm)
}

.sc_aux_matches <- function(files, prefix, kind) {
  nm <- tolower(basename(files))
  nm <- sub("\\.gz$", "", nm)
  suffix <- switch(kind,
    features = c("features.tsv", "genes.tsv", "features.csv", "genes.csv"),
    barcodes = c("barcodes.tsv", "barcodes.csv")
  )
  hits <- files[vapply(nm, function(x) {
    any(x == suffix) || (nzchar(prefix) && any(x == paste0(prefix, "_", suffix)))
  }, logical(1))]
  unique(hits)
}

.sc_scan_table_hint <- function(path) {
  lines <- tryCatch(.sc_read_lines(path, n = 2L), error = function(e) character())
  if (!length(lines)) return(list(sep = NA_character_, orientation = "features_by_cells", confidence = "default", reason = "could not inspect header", matrix_like = FALSE))
  sep <- .sc_detect_separator(lines[[1L]])
  fields <- strsplit(lines[[1L]], sep, fixed = TRUE)[[1L]]
  second <- if (length(lines) >= 2L) strsplit(lines[[2L]], sep, fixed = TRUE)[[1L]] else character()
  if (length(fields) < 3L || length(second) < 3L) {
    return(list(sep = sep, orientation = "features_by_cells", confidence = "default", reason = "header does not look like a matrix", matrix_like = FALSE))
  }
  numeric_values <- suppressWarnings(as.numeric(second[-1L]))
  if (mean(!is.na(numeric_values)) < 0.8) {
    return(list(sep = sep, orientation = "features_by_cells", confidence = "default", reason = "expression columns are not predominantly numeric", matrix_like = FALSE))
  }
  hint <- .sc_guess_orientation(second[[1L]], fields[-1L], fields[[1L]])
  list(sep = sep, orientation = hint$value, confidence = hint$confidence, reason = hint$reason, matrix_like = TRUE)
}

.sc_scan_one <- function(path, root) {
  fmt <- .sc_detect_table_format(path)
  if (is.na(fmt)) return(NULL)
  nm <- tolower(basename(path))
  if (fmt == "10x_mtx_component") return(NULL)
  if (fmt %in% c("csv", "tsv", "txt") &&
      grepl("(^|_)(features|genes|barcodes|metadata|annotation|sample[_-]?sheet)(_|\\.)", nm)) return(NULL)

  sample_id <- .sc_sample_id(path, root)
  hint <- list(sep = NA_character_, orientation = "features_by_cells", confidence = "not_applicable", reason = "format defines orientation")
  status <- "candidate"
  note <- ""
  if (fmt %in% c("csv", "tsv", "txt")) {
    hint <- .sc_scan_table_hint(path)
    if (!hint$matrix_like) return(NULL)
    if (hint$confidence == "default") note <- hint$reason
  } else if (fmt == "hdf5") {
    fmt <- "10x_h5"
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      status <- "needs_seurat"
      note <- "Install Seurat to read 10x HDF5 files."
    }
  } else if (fmt == "h5ad") {
    status <- "unsupported_yet"
    note <- "H5AD import is planned for a later release."
  } else if (fmt == "rds") {
    status <- "class_checked_on_read"
  }
  data.frame(
    dataset_id = sample_id,
    sample_id = sample_id,
    format = fmt,
    matrix_path = normalizePath(path, winslash = "/", mustWork = FALSE),
    features_path = NA_character_,
    barcodes_path = NA_character_,
    orientation = hint$orientation,
    orientation_confidence = hint$confidence,
    orientation_reason = hint$reason,
    status = status,
    note = note,
    stringsAsFactors = FALSE
  )
}

#' Scan a directory for single-cell datasets
#'
#' Inspects file names, gzip signatures, table headers, and 10x sidecar files.
#' Scanning does not read complete expression matrices.
#'
#' @param path A file or directory to inspect.
#' @param recursive Whether to scan nested directories.
#' @param orientation Logical fallback orientation: `TRUE` means features in
#'   rows and cells in columns; `FALSE` means cells in rows and features in columns.
#' @return A data frame with one row per candidate dataset and class `sc_scan`.
#' @export
scan_sc_data <- function(path, recursive = TRUE,
                         orientation = TRUE) {
  orientation <- .sc_normalize_orientation(orientation)
  if (length(path) != 1L) .sc_abort("`path` must be one file or directory.")
  if (!length(path) || is.na(path) || !file.exists(path)) .sc_abort("Path does not exist: ", path)
  selected_file <- NULL
  if (dir.exists(path)) {
    root <- normalizePath(path, winslash = "/", mustWork = TRUE)
    files <- list.files(root, recursive = recursive, full.names = TRUE,
                        all.files = FALSE, include.dirs = FALSE, no.. = TRUE)
  } else {
    root <- dirname(normalizePath(path, winslash = "/", mustWork = TRUE))
    selected_file <- normalizePath(path, winslash = "/", mustWork = TRUE)
    if (grepl("\\.mtx(\\.gz)?$", tolower(selected_file))) {
      files <- list.files(root, full.names = TRUE, all.files = FALSE,
                          include.dirs = FALSE, no.. = TRUE)
      files <- normalizePath(files, winslash = "/", mustWork = FALSE)
    } else {
      files <- selected_file
    }
  }
  files <- files[file.info(files)$isdir %in% FALSE]
  rows <- list()

  # Group Matrix Market files with their feature and barcode sidecars.
  mtx <- files[grepl("\\.mtx(\\.gz)?$", tolower(files))]
  for (mp in mtx) {
    stem <- .sc_mtx_stem(mp)
    prefix <- sub("(_matrix|matrix)$", "", stem, ignore.case = TRUE)
    siblings <- files[dirname(files) == dirname(mp)]
    f <- .sc_aux_matches(siblings, prefix, "features")
    b <- .sc_aux_matches(siblings, prefix, "barcodes")
    status <- if (length(f) == 1L && length(b) == 1L) "candidate" else "needs_sidecars"
    note <- if (status == "candidate") "" else paste0(
      "Expected one features/genes file and one barcodes file; found ",
      length(f), " and ", length(b), "."
    )
    sid <- .sc_sample_id(mp, root)
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_id = sid, sample_id = sid, format = "10x_mtx",
      matrix_path = normalizePath(mp, winslash = "/", mustWork = FALSE),
      features_path = if (length(f) == 1L) normalizePath(f, winslash = "/", mustWork = FALSE) else NA_character_,
      barcodes_path = if (length(b) == 1L) normalizePath(b, winslash = "/", mustWork = FALSE) else NA_character_,
      orientation = "features_by_cells", orientation_confidence = "high",
      orientation_reason = "10x Matrix Market convention", status = status,
      note = note, stringsAsFactors = FALSE
    )
  }

  other <- files[!grepl("\\.mtx(\\.gz)?$", tolower(files))]
  for (f in other) {
    row <- .sc_scan_one(f, root)
    if (!is.null(row)) rows[[length(rows) + 1L]] <- row
  }

  if (!length(rows)) {
    out <- data.frame(
      dataset_id = character(), sample_id = character(), format = character(),
      matrix_path = character(), features_path = character(), barcodes_path = character(),
      orientation = character(), orientation_confidence = character(),
      orientation_reason = character(), status = character(), note = character(),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, rows)
    # A generic text file that is part of a discovered 10x triplet is not a
    # second independent expression matrix.
    sidecars <- unique(c(out$features_path, out$barcodes_path))
    sidecars <- sidecars[!is.na(sidecars)]
    out <- out[!out$matrix_path %in% sidecars, , drop = FALSE]
    defaulted <- out$orientation_confidence == "default"
    out$orientation[defaulted] <- .sc_orientation_label(orientation)
    out$orientation_reason[defaulted] <- paste0("ambiguous; default orientation = ", orientation)
    if (!is.null(selected_file)) out <- out[out$matrix_path == selected_file, , drop = FALSE]
    rownames(out) <- NULL
  }
  attr(out, "root") <- root
  class(out) <- c("sc_scan", class(out))
  out
}

#' @export
print.sc_scan <- function(x, ...) {
  cat("scAutoLoad scan: ", nrow(x), " candidate dataset(s)\n", sep = "")
  if (!nrow(x)) {
    cat("No supported files were found.\n")
    return(invisible(x))
  }
  show <- x[, c("sample_id", "format", "orientation", "orientation_confidence", "status", "note"), drop = FALSE]
  # Subsetting a data frame preserves the `sc_scan` class, which would dispatch
  # back to this method and recurse indefinitely. Call the data-frame method.
  print.data.frame(show, row.names = FALSE, ...)
  invisible(x)
}

