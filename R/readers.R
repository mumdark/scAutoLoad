.sc_read_aux_table <- function(path) {
  first <- .sc_read_lines(path, n = 1L)
  sep <- .sc_detect_separator(if (length(first)) first[[1L]] else "")
  con <- .sc_open_text(path)
  on.exit(close(con), add = TRUE)
  utils::read.table(con, header = FALSE, sep = sep, quote = "", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE,
                    colClasses = "character")
}

.sc_read_mtx <- function(candidate) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    .sc_abort("Package 'Matrix' is required to read Matrix Market files.")
  }
  if (is.na(candidate$features_path) || is.na(candidate$barcodes_path)) {
    .sc_abort("This 10x MTX dataset is missing its features/genes or barcodes sidecar. See `scan_sc_data()` for details.")
  }
  con <- .sc_open_text(candidate$matrix_path)
  mat <- tryCatch(Matrix::readMM(con), finally = close(con))
  feature_tab <- .sc_read_aux_table(candidate$features_path)
  barcode_tab <- .sc_read_aux_table(candidate$barcodes_path)
  if (nrow(feature_tab) != nrow(mat)) {
    .sc_abort("Feature file has ", nrow(feature_tab), " rows but MTX has ", nrow(mat), " rows.")
  }
  if (nrow(barcode_tab) != ncol(mat)) {
    .sc_abort("Barcode file has ", nrow(barcode_tab), " rows but MTX has ", ncol(mat), " columns.")
  }
  feature_col <- if (ncol(feature_tab) >= 2L) 2L else 1L
  feature_names <- .sc_make_unique_names(feature_tab[[feature_col]], "feature")
  cell_names <- .sc_make_unique_names(barcode_tab[[1L]], "cell")
  rownames(mat) <- feature_names$names
  colnames(mat) <- cell_names$names
  list(matrix = mat, orientation = "features_by_cells", orientation_confidence = "high",
       orientation_reason = "10x Matrix Market convention", data_type = "counts",
       warnings = character(), original_feature_names = feature_names$original,
       original_cell_names = cell_names$original)
}

.sc_read_table <- function(path, delimiter = NULL, orientation = TRUE,
                           auto_orientation = TRUE, na_to_zero = TRUE) {
  lines <- .sc_read_lines(path, n = 1L)
  sep <- .sc_detect_separator(if (length(lines)) lines[[1L]] else "", explicit = delimiter)
  gzip <- .sc_is_gzip(path)
  suffix_gzip <- grepl("\\.gz$", tolower(path))
  use_fread <- requireNamespace("data.table", quietly = TRUE) && identical(gzip, suffix_gzip)
  if (use_fread) {
    dat <- data.table::fread(
      file = path, sep = sep, header = TRUE, data.table = FALSE,
      check.names = FALSE, na.strings = "NA", showProgress = interactive()
    )
  } else {
    con <- .sc_open_text(path)
    on.exit(close(con), add = TRUE)
    dat <- utils::read.table(
      con, header = TRUE, row.names = NULL, sep = sep, quote = "\"", comment.char = "",
      check.names = FALSE, stringsAsFactors = FALSE, fill = FALSE,
      blank.lines.skip = TRUE, na.strings = "NA"
    )
  }
  if (ncol(dat) < 2L) .sc_abort("The table needs a feature/cell identifier column and at least one expression column.")
  first_col <- dat[[1L]]
  ids <- as.character(first_col)
  cell_or_feature_names <- names(dat)[-1L]
  raw <- dat[-1L]
  missing_values_replaced <- 0L
  numeric_cols <- lapply(seq_along(raw), function(i) {
    original <- raw[[i]]
    converted <- suppressWarnings(as.numeric(as.character(original)))
    bad <- is.na(converted) & !is.na(original) & nzchar(trimws(as.character(original)))
    if (any(bad)) {
      .sc_abort("Non-numeric expression value in column '", names(raw)[i],
                "'. Check whether this is an expression matrix with a single identifier column.")
    }
    if (isTRUE(na_to_zero)) {
      missing <- is.na(converted)
      if (any(missing)) {
        missing_values_replaced <<- missing_values_replaced + sum(missing)
        converted[missing] <- 0
      }
    }
    converted
  })
  mat <- do.call(cbind, numeric_cols)
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = length(numeric_cols))

  inferred <- .sc_guess_orientation(ids, cell_or_feature_names,
                                    row_header = names(dat)[[1L]])
  if (isTRUE(auto_orientation) && inferred$confidence != "default") {
    resolved <- inferred$value
    confidence <- inferred$confidence
    reason <- inferred$reason
  } else if (isTRUE(auto_orientation)) {
    resolved <- .sc_orientation_label(orientation)
    confidence <- "default"
    reason <- paste0("ambiguous; used orientation = ", orientation)
  } else {
    resolved <- .sc_orientation_label(orientation)
    confidence <- "user"
    reason <- "user-specified orientation"
  }
  if (resolved == "cells_by_features") mat <- t(mat)
  row_names <- if (resolved == "features_by_cells") {
    .sc_make_unique_names(ids, "feature")
  } else {
    .sc_make_unique_names(cell_or_feature_names, "feature")
  }
  col_names <- if (resolved == "features_by_cells") {
    .sc_make_unique_names(cell_or_feature_names, "cell")
  } else {
    .sc_make_unique_names(ids, "cell")
  }
  rownames(mat) <- row_names$names
  colnames(mat) <- col_names$names
  original_features <- row_names$original
  original_cells <- col_names$original

  warnings <- character()
  if (isTRUE(auto_orientation) && identical(inferred$confidence, "default")) {
    warnings <- c(warnings, paste0("Orientation was ambiguous; used orientation = ", orientation, "."))
  }
  list(matrix = mat, orientation = resolved, orientation_confidence = confidence,
       orientation_reason = reason, data_type = "counts", warnings = warnings,
       missing_values_replaced = missing_values_replaced,
       original_feature_names = original_features, original_cell_names = original_cells,
       delimiter = sep)
}

.sc_read_h5 <- function(path) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    .sc_abort("Reading 10x HDF5 requires the optional 'Seurat' package.")
  }
  x <- Seurat::Read10X_h5(path)
  warnings <- character()
  if (is.list(x)) {
    nm <- names(x)
    idx <- match("Gene Expression", nm)
    if (is.na(idx)) {
      .sc_abort("The 10x HDF5 file does not contain a 'Gene Expression' matrix. Available modalities: ",
                paste(nm, collapse = ", "))
    }
    if (length(x) > 1L) {
      warnings <- paste0("Additional 10x modalities were not imported in this release: ",
                         paste(setdiff(nm, "Gene Expression"), collapse = ", "))
    }
    x <- x[[idx]]
  }
  list(matrix = x, orientation = "features_by_cells", orientation_confidence = "high",
       orientation_reason = "10x HDF5 convention", data_type = "counts",
       warnings = warnings, original_feature_names = rownames(x))
}

.sc_matrix_from_rds <- function(path) {
  x <- readRDS(path)
  if (inherits(x, "Seurat")) {
    return(list(object = x, source_class = class(x)))
  }
  if (inherits(x, "SingleCellExperiment")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      .sc_abort("Converting a SingleCellExperiment RDS requires 'SeuratObject'.")
    }
    return(list(object = SeuratObject::as.Seurat(x), source_class = class(x)))
  }
  if (is.matrix(x) || inherits(x, "Matrix")) {
    return(list(matrix = x, orientation = "features_by_cells",
                orientation_confidence = "assumed", orientation_reason = "R matrix dimnames convention",
                data_type = "counts", warnings = character(),
                original_feature_names = rownames(x), original_cell_names = colnames(x),
                source_class = class(x)))
  }
  .sc_abort("Unsupported RDS class: ", paste(class(x), collapse = "/"),
            ". Supported RDS inputs are Seurat, SingleCellExperiment, and numeric matrices.")
}

.sc_read_matrix_source <- function(candidate, orientation, delimiter,
                                   auto_orientation = TRUE, na_to_zero = TRUE) {
  switch(candidate$format,
    `10x_mtx` = .sc_read_mtx(candidate),
    `10x_h5` = .sc_read_h5(candidate$matrix_path),
    h5ad = .sc_abort("H5AD import is planned for a later release."),
    csv = .sc_read_table(candidate$matrix_path, delimiter, orientation, auto_orientation, na_to_zero),
    tsv = .sc_read_table(candidate$matrix_path, delimiter %||% "\t", orientation, auto_orientation, na_to_zero),
    txt = .sc_read_table(candidate$matrix_path, delimiter, orientation, auto_orientation, na_to_zero),
    rds = .sc_matrix_from_rds(candidate$matrix_path),
    .sc_abort("Unsupported input format: ", candidate$format)
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

