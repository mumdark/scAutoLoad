.sc_abort <- function(..., call. = FALSE) {
  stop(..., call. = call.)
}

.sc_is_gzip <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 2L)
  length(magic) == 2L && identical(magic, as.raw(c(31L, 139L)))
}

.sc_open_text <- function(path) {
  if (.sc_is_gzip(path)) {
    gzfile(path, open = "rt", encoding = "UTF-8")
  } else {
    file(path, open = "rt", encoding = "UTF-8")
  }
}

.sc_read_lines <- function(path, n = -1L) {
  con <- .sc_open_text(path)
  on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE, encoding = "UTF-8")
}

.sc_detect_table_format <- function(path) {
  nm <- tolower(basename(path))
  if (grepl("\\.mtx(\\.gz)?$", nm)) return("10x_mtx_component")
  if (grepl("\\.h5ad$", nm)) return("h5ad")
  if (grepl("\\.h5$|\\.hdf5$", nm)) return("hdf5")
  if (grepl("\\.rds$", nm)) return("rds")
  if (grepl("\\.(csv|tsv|txt)(\\.gz)?$", nm) ||
      grepl("\\.gz\\.(csv|tsv|txt)$", nm)) {
    if (grepl("\\.csv(\\.gz)?$|\\.gz\\.csv$", nm)) return("csv")
    if (grepl("\\.tsv(\\.gz)?$|\\.gz\\.tsv$", nm)) return("tsv")
    return("txt")
  }
  NA_character_
}

.sc_detect_separator <- function(line, explicit = NULL) {
  if (!is.null(explicit)) {
    if (!explicit %in% c(",", "\t", ";")) {
      .sc_abort("`delimiter` must be one of ',', '\\t', or ';'.")
    }
    return(explicit)
  }
  if (!length(line) || is.na(line) || !nzchar(line)) return("\t")
  candidates <- c("\t", ",", ";")
  counts <- vapply(candidates, function(sep) {
    lengths(regmatches(line, gregexpr(sep, line, fixed = TRUE)))
  }, integer(1))
  if (max(counts) == 0L) return("\t")
  candidates[[which.max(counts)]]
}

.sc_feature_score <- function(x) {
  x <- tolower(as.character(x))
  if (!length(x)) return(0)
  mean(grepl("gene|feature|symbol|ensembl|^id$", x)) * 3
}

.sc_cell_score <- function(x) {
  x <- tolower(as.character(x))
  if (!length(x)) return(0)
  mean(grepl("cell|barcode|^bc$|cell_id|barcode_id", x)) * 3
}

.sc_barcode_fraction <- function(x) {
  x <- as.character(x)
  if (!length(x)) return(0)
  mean(grepl("^[ACGTN]{10,}(-[0-9]+)?$", x, ignore.case = TRUE))
}

.sc_feature_id_fraction <- function(x) {
  x <- as.character(x)
  if (!length(x)) return(0)
  mean(grepl("^(ENSG|ENSMUSG|ENSDARG|ENSGALG)[0-9]+", x, ignore.case = TRUE))
}

.sc_guess_orientation <- function(row_ids, col_ids, row_header = "") {
  row_ids <- head(as.character(row_ids), 200L)
  col_ids <- head(as.character(col_ids), 200L)
  row_header <- tolower(as.character(row_header)[1L])
  score_features_in_rows <- .sc_feature_score(row_header) +
    3 * .sc_feature_id_fraction(row_ids) - 3 * .sc_barcode_fraction(row_ids) +
    3 * .sc_barcode_fraction(col_ids) + .sc_cell_score(col_ids)
  score_features_in_cols <- .sc_cell_score(row_header) +
    3 * .sc_barcode_fraction(row_ids) + .sc_cell_score(row_ids) +
    3 * .sc_feature_id_fraction(col_ids) - 3 * .sc_barcode_fraction(col_ids)

  delta <- score_features_in_rows - score_features_in_cols
  if (delta >= 2) {
    list(value = "features_by_cells", confidence = "high", reason = "feature-like row identifiers or cell-like column identifiers")
  } else if (delta <= -2) {
    list(value = "cells_by_features", confidence = "high", reason = "cell-like row identifiers or feature-like column identifiers")
  } else if (delta != 0) {
    list(value = if (delta > 0) "features_by_cells" else "cells_by_features", confidence = "medium", reason = "weak identifier evidence")
  } else {
    list(value = "features_by_cells", confidence = "default", reason = "ambiguous; used the default features-in-rows orientation")
  }
}

.sc_sample_id <- function(path, root = dirname(path)) {
  parent <- basename(dirname(path))
  if (!identical(normalizePath(dirname(path), winslash = "/", mustWork = FALSE),
                 normalizePath(root, winslash = "/", mustWork = FALSE))) {
    if (!parent %in% c("filtered_feature_bc_matrix", "raw_feature_bc_matrix")) return(parent)
    return(basename(dirname(dirname(path))))
  }
  nm <- basename(path)
  nm <- sub("\\.gz$", "", nm, ignore.case = TRUE)
  nm <- sub("\\.(csv|tsv|txt|mtx|h5|hdf5|rds)$", "", nm, ignore.case = TRUE)
  nm <- sub("(_filtered_feature_bc_matrix|_raw_feature_bc_matrix|_matrix|_counts|_count|_expression|_data)$", "", nm, ignore.case = TRUE)
  nm <- sub("^(filtered_feature_bc_matrix|raw_feature_bc_matrix|matrix|counts|count|expression|data)$", basename(root), nm, ignore.case = TRUE)
  if (!nzchar(nm)) basename(root) else nm
}

.sc_make_unique_names <- function(x, what = "names") {
  x <- as.character(x)
  blank <- is.na(x) | !nzchar(x)
  if (any(blank)) x[blank] <- paste0(what, "_", which(blank))
  if (anyDuplicated(x)) {
    original <- x
    x <- make.unique(x, sep = "_")
    return(list(names = x, original = original))
  }
  list(names = x, original = x)
}

.sc_to_sparse <- function(x) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    .sc_abort("Package 'Matrix' is required to construct a sparse expression matrix.")
  }
  if (inherits(x, "sparseMatrix")) return(x)
  Matrix::Matrix(x, sparse = TRUE)
}

.sc_value_flags <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    vals <- tryCatch(x@x, error = function(e) rep(1, nrow(Matrix::summary(x))))
    return(c(
      finite = all(is.finite(vals)),
      nonnegative = all(vals >= 0, na.rm = TRUE),
      integerish = all(abs(vals - round(vals)) <= 1e-8, na.rm = TRUE)
    ))
  }
  finite <- nonnegative <- integerish <- TRUE
  for (j in seq_len(ncol(x))) {
    vals <- x[, j]
    finite <- finite && all(is.finite(vals))
    nonnegative <- nonnegative && all(vals >= 0, na.rm = TRUE)
    integerish <- integerish && all(abs(vals - round(vals)) <= 1e-8, na.rm = TRUE)
  }
  c(finite = finite, nonnegative = nonnegative, integerish = integerish)
}

