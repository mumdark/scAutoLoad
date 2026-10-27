.sc_normalize_orientation <- function(orientation) {
  if (!is.logical(orientation) || length(orientation) != 1L || is.na(orientation)) {
    .sc_abort("`orientation` must be one logical value: TRUE (features in rows) or FALSE (features in columns).")
  }
  orientation
}

.sc_orientation_label <- function(orientation) {
  if (isTRUE(orientation)) "features_by_cells" else "cells_by_features"
}

.sc_candidate_from_path <- function(path, format = "auto", recursive = TRUE) {
  found <- scan_sc_data(path, recursive = recursive)
  if (!nrow(found)) .sc_abort("No supported single-cell dataset was found under: ", path)
  if (nrow(found) > 1L) {
    .sc_abort("Found ", nrow(found), " candidate datasets. Use `read_sc_batch()` or pass a path containing one dataset.")
  }
  if (found$status[[1L]] == "needs_sidecars") {
    .sc_abort(found$note[[1L]], " Run `scan_sc_data()` to inspect the missing 10x files.")
  }
  if (found$status[[1L]] == "unsupported_yet") {
    .sc_abort(found$note[[1L]])
  }
  if (!identical(format, "auto")) found$format[[1L]] <- format
  as.list(found[1L, , drop = FALSE])
}

.sc_normalize_import_names <- function(imported) {
  mat <- imported$matrix
  if (is.null(dim(mat)) || length(dim(mat)) != 2L || any(dim(mat) == 0L)) {
    .sc_abort("Expression input must be a non-empty two-dimensional matrix.")
  }
  original_features <- imported$original_feature_names %||% rownames(mat)
  original_cells <- imported$original_cell_names %||% colnames(mat)
  if (is.null(original_features)) original_features <- paste0("feature_", seq_len(nrow(mat)))
  if (is.null(original_cells)) original_cells <- paste0("cell_", seq_len(ncol(mat)))
  if (length(original_features) != nrow(mat) || length(original_cells) != ncol(mat)) {
    .sc_abort("Feature/cell labels do not match the expression matrix dimensions.")
  }
  feature_names <- .sc_make_unique_names(original_features, "feature")
  cell_names <- .sc_make_unique_names(original_cells, "cell")
  rownames(mat) <- feature_names$names
  colnames(mat) <- cell_names$names
  imported$matrix <- mat
  imported$original_feature_names <- feature_names$original
  imported$original_cell_names <- cell_names$original
  imported$duplicate_features_repaired <- anyDuplicated(original_features) > 0L
  imported$duplicate_cells_repaired <- anyDuplicated(original_cells) > 0L
  imported
}

.sc_check_counts <- function(mat, data_type) {
  warnings <- character()
  flags <- .sc_value_flags(mat)
  if (!flags[["finite"]]) {
    .sc_abort("The matrix contains NA, NaN, or infinite values. Clean the expression matrix before import.")
  }
  if (data_type == "counts") {
    if (!flags[["nonnegative"]]) {
      .sc_abort("The matrix contains negative values and cannot be imported as counts. Use `data_type = 'data'` for normalized or scaled expression.")
    }
    if (!flags[["integerish"]]) {
      warnings <- c(warnings, "The selected Counts matrix contains non-integer values; confirm that it contains raw counts.")
    }
  }
  warnings
}

.sc_replace_missing_values <- function(mat) {
  if (inherits(mat, "sparseMatrix") && "x" %in% methods::slotNames(mat)) {
    values <- methods::slot(mat, "x")
    missing <- is.na(values)
    n_missing <- sum(missing)
    if (n_missing) {
      values[missing] <- 0
      methods::slot(mat, "x") <- values
      mat <- Matrix::drop0(mat)
    }
    return(list(matrix = mat, n_missing = n_missing))
  }
  missing <- is.na(mat)
  n_missing <- sum(missing)
  if (n_missing) mat[missing] <- 0
  list(matrix = mat, n_missing = n_missing)
}

#' Read a single-cell dataset into a Seurat object
#'
#' @param path A matrix file, RDS/HDF5 file, or directory containing one dataset.
#' @param sample_id Optional sample identifier; defaults to a name inferred
#'   from the file or containing folder and is stored in cell metadata.
#' @param format Input format, or `"auto"` to detect it.
#' @param version Seurat assay layout: `"auto"`, `"v4"` (legacy Assay), or `"v5"` (Assay5).
#' @param orientation Logical orientation: `TRUE` means features in rows and
#'   cells in columns; `FALSE` means cells in rows and features in columns. If
#'   omitted, the function detects orientation from identifiers and falls back
#'   to `TRUE` when ambiguous. Explicitly passing `TRUE` or `FALSE` forces it.
#' @param data_type `"counts"` by default; use `"data"` for normalized expression.
#' @param na_to_zero If `TRUE` (default), replace missing expression values with
#'   zero and print how many values were replaced. Set `FALSE` to keep them and
#'   stop with an error before creating the Seurat object.
#' @param assay Assay name for the expression matrix.
#' @param project Project name stored in the Seurat object.
#' @param metadata Optional cell-level data frame or delimited file. Cell names
#'   must be row names or provided in a `cell` column; a single row is treated
#'   as sample-level metadata and repeated for all cells.
#' @param delimiter Optional delimiter override for text matrices: comma, tab, or semicolon.
#' @param recursive Scan nested directories when `path` is a directory.
#' @param return_report If `TRUE`, return `list(object, report)`.
#' @return A Seurat object, or an object and import report when requested.
#' @export
read_sc_data <- function(path, sample_id = NULL, format = "auto",
                         version = "auto",
                         orientation = TRUE,
                         data_type = "counts",
                         assay = "RNA", project = NULL,
                         metadata = NULL, delimiter = NULL,
                         recursive = TRUE, return_report = FALSE,
                         na_to_zero = TRUE) {
  orientation_was_missing <- missing(orientation)
  orientation <- .sc_normalize_orientation(orientation)
  if (!is.logical(na_to_zero) || length(na_to_zero) != 1L || is.na(na_to_zero)) {
    .sc_abort("`na_to_zero` must be one logical value: TRUE or FALSE.")
  }
  data_type <- match.arg(data_type, c("counts", "data"))
  version <- match.arg(version, c("auto", "v4", "v5"))
  if (length(path) != 1L || !file.exists(path)) .sc_abort("`path` must point to an existing file or directory.")
  candidate <- .sc_candidate_from_path(path, format = format, recursive = recursive)
  fmt <- candidate$format[[1L]]
  inferred_sample <- candidate$sample_id[[1L]]
  if (!is.null(sample_id)) inferred_sample <- as.character(sample_id)[1L]
  if (is.null(project)) project <- inferred_sample

  if (fmt %in% c("csv", "tsv", "txt")) {
    source_size <- file.info(candidate$matrix_path[[1L]])$size
    if (is.finite(source_size) && source_size >= 50 * 1024^2) {
      warning("This dense text matrix is at least 50 MiB and may require several times its file size in memory during import.", call. = FALSE)
    }
  }

  imported <- .sc_read_matrix_source(candidate, orientation, delimiter,
                                     auto_orientation = orientation_was_missing,
                                     na_to_zero = na_to_zero)
  if (!is.null(imported$object)) {
    object <- imported$object
    if (!is.null(metadata)) {
      if (!requireNamespace("SeuratObject", quietly = TRUE)) .sc_require_seuratobject()
      metadata_df <- .sc_read_metadata(metadata, colnames(object), inferred_sample)
      object <- SeuratObject::AddMetaData(object, metadata = metadata_df)
    }
    object <- .sc_set_sample_identity(object, inferred_sample)
    report <- list(
      sample_id = inferred_sample, source = candidate$matrix_path[[1L]],
      format = fmt, source_class = imported$source_class,
      orientation = "preserved_from_object", data_type = "preserved_from_object",
      assay_version = paste(class(object[[SeuratObject::DefaultAssay(object)]]), collapse = "/"),
      warnings = character()
    )
    object@misc$scAutoLoad <- report
    if (return_report) return(list(object = object, report = report))
    return(object)
  }

  if (is.null(imported$missing_values_replaced)) {
    imported$missing_values_replaced <- 0L
    if (isTRUE(na_to_zero)) {
      missing_result <- .sc_replace_missing_values(imported$matrix)
      imported$matrix <- missing_result$matrix
      imported$missing_values_replaced <- missing_result$n_missing
    }
  }
  if (isTRUE(na_to_zero) && imported$missing_values_replaced > 0L) {
    message("scAutoLoad: replaced ", format(imported$missing_values_replaced, big.mark = ","),
            " missing expression value(s) with 0 (na_to_zero = TRUE).")
  }

  if (fmt %in% c("10x_mtx", "10x_h5") && !orientation_was_missing) {
    resolved <- .sc_orientation_label(orientation)
    if (!orientation) imported$matrix <- t(imported$matrix)
    imported$orientation <- resolved
    imported$orientation_confidence <- "user"
    imported$orientation_reason <- "user-specified orientation"
  } else if (fmt == "rds" && orientation_was_missing) {
    rowids <- rownames(imported$matrix)
    colids <- colnames(imported$matrix)
    guess <- .sc_guess_orientation(rowids, colids, "")
    if (guess$value == "cells_by_features") imported$matrix <- t(imported$matrix)
    imported$orientation <- guess$value
    imported$orientation_confidence <- guess$confidence
    imported$orientation_reason <- guess$reason
    if (identical(guess$confidence, "default")) {
      imported$orientation <- .sc_orientation_label(orientation)
      imported$warnings <- c(imported$warnings,
                             "Orientation was ambiguous; used default TRUE (features in rows, cells in columns).")
    }
  } else if (fmt == "rds" && !orientation_was_missing) {
    resolved <- .sc_orientation_label(orientation)
    if (!orientation) imported$matrix <- t(imported$matrix)
    imported$orientation <- resolved
    imported$orientation_confidence <- "user"
    imported$orientation_reason <- "user-specified orientation"
  }

  imported <- .sc_normalize_import_names(imported)
  warnings <- unique(c(imported$warnings, .sc_check_counts(imported$matrix, data_type)))
  if (imported$duplicate_features_repaired) {
    warnings <- c(warnings, "Duplicate feature names were made unique; duplicate rows were not summed.")
  }
  if (imported$duplicate_cells_repaired) {
    warnings <- c(warnings, "Duplicate cell names within this dataset were made unique.")
  }

  message(
    "scAutoLoad: orientation = ", imported$orientation,
    " (", if (imported$orientation == "features_by_cells") "features in rows, cells in columns" else "cells in rows, features in columns",
    "; confidence: ", imported$orientation_confidence, "; ", imported$orientation_reason, ")."
  )
  if (length(warnings)) for (w in warnings) warning(w, call. = FALSE)

  metadata_df <- .sc_read_metadata(metadata, colnames(imported$matrix), inferred_sample)
  created <- .sc_create_seurat(
    imported$matrix, data_type = data_type, assay = assay,
    project = project, metadata = metadata_df, version = version
  )
  report <- list(
    sample_id = inferred_sample,
    source = candidate$matrix_path[[1L]],
    format = fmt,
    dimensions_features_cells = dim(imported$matrix),
    orientation = imported$orientation,
    orientation_confidence = imported$orientation_confidence,
    orientation_reason = imported$orientation_reason,
    data_type = data_type,
    missing_values_replaced = imported$missing_values_replaced,
    delimiter = imported$delimiter %||% NA_character_,
    assay = assay,
    assay_version = created$assay_version,
    seuratobject_version = created$seuratobject_version,
    original_feature_names = imported$original_feature_names,
    original_cell_names = imported$original_cell_names,
    warnings = warnings
  )
  object <- .sc_set_sample_identity(created$object, inferred_sample)
  object@misc$scAutoLoad <- report
  if (return_report) return(list(object = object, report = report))
  object
}
