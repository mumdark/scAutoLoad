.sc_sample_metadata_row <- function(sample_metadata, sample_id) {
  if (is.null(sample_metadata)) return(NULL)
  if (!is.data.frame(sample_metadata)) .sc_abort("`sample_metadata` must be a data frame.")
  if ("sample_id" %in% names(sample_metadata)) {
    idx <- which(as.character(sample_metadata$sample_id) == sample_id)
    if (length(idx) != 1L) .sc_abort("`sample_metadata` must contain exactly one row for sample '", sample_id, "'.")
    out <- sample_metadata[idx, setdiff(names(sample_metadata), "sample_id"), drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  if (!is.null(rownames(sample_metadata)) && sample_id %in% rownames(sample_metadata)) {
    out <- sample_metadata[sample_id, , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  .sc_abort("`sample_metadata` needs a 'sample_id' column or sample IDs as row names.")
}

#' Read all detected samples in a directory
#'
#' @param path Directory to scan.
#' @param recursive Whether to scan nested directories.
#' @param sample_metadata Optional sample-level data frame keyed by `sample_id`.
#' @param merge Combine samples into one Seurat object (default). Sample names
#'   inferred from file or folder names are prefixed to cell barcodes and stored
#'   in the `sample_id` and `orig.ident` metadata columns. Set `FALSE` to return
#'   a named list of per-sample objects.
#' @param version Seurat assay layout, passed to [read_sc_data()].
#' @param orientation Logical orientation: `TRUE` means features in rows and
#'   cells in columns; `FALSE` means cells in rows and features in columns. If
#'   omitted, each sample is auto-detected and ambiguous matrices use `TRUE`.
#' @param data_type `"counts"` by default; use `"data"` for normalized expression.
#' @param na_to_zero If `TRUE` (default), replace missing expression values with
#'   zero and print the number replaced for each sample. Set `FALSE` to retain
#'   them and stop before creating the Seurat object.
#' @param assay Assay name.
#' @param delimiter Optional delimiter override for text files.
#' @param return_report If `TRUE`, return imported object(s) with a combined report.
#' @param formats Optional character vector selecting formats from the scan,
#'   such as `"10x_mtx"`, `"csv"`, `"tsv"`, `"txt"`, `"10x_h5"`, or `"rds"`.
#' @return By default, a merged Seurat object. With `merge = FALSE`, a named
#'   list of per-sample objects. With `return_report = TRUE`, the object(s) and
#'   import reports are returned in a list.
#' @export
read_sc_batch <- function(path, recursive = TRUE, sample_metadata = NULL,
                          merge = TRUE, version = "auto",
                          orientation = TRUE, data_type = "counts",
                          assay = "RNA", delimiter = NULL,
                          return_report = FALSE, formats = NULL,
                          na_to_zero = TRUE) {
  orientation_was_missing <- missing(orientation)
  orientation <- .sc_normalize_orientation(orientation)
  if (!is.logical(na_to_zero) || length(na_to_zero) != 1L || is.na(na_to_zero)) {
    .sc_abort("`na_to_zero` must be one logical value: TRUE or FALSE.")
  }
  if (!is.logical(merge) || length(merge) != 1L || is.na(merge)) {
    .sc_abort("`merge` must be one logical value: TRUE or FALSE.")
  }
  found <- scan_sc_data(path, recursive = recursive, orientation = orientation)
  if (!is.null(formats)) {
    supported <- c("10x_mtx", "csv", "tsv", "txt", "10x_h5", "rds")
    if (!is.character(formats) || !length(formats) || anyNA(formats) || any(!formats %in% supported)) {
      .sc_abort("`formats` must contain supported formats: ", paste(supported, collapse = ", "), ".")
    }
    found <- found[found$format %in% formats, , drop = FALSE]
  }
  if (!nrow(found)) .sc_abort("No supported datasets were found under: ", path)
  incomplete <- found$status == "needs_sidecars"
  if (any(incomplete)) {
    .sc_abort("One or more 10x MTX datasets have missing sidecars: ",
              paste(found$sample_id[incomplete], collapse = ", "),
              ". Resolve these first or remove them from the input directory.")
  }
  ids <- as.character(found$sample_id)
  if (anyDuplicated(ids)) {
    .sc_abort("Sample IDs are not unique: ", paste(unique(ids[duplicated(ids)]), collapse = ", "),
              ". Put samples in separate folders or rename files so each sample has a unique ID.")
  }
  objects <- vector("list", nrow(found))
  reports <- vector("list", nrow(found))
  names(objects) <- ids
  for (i in seq_len(nrow(found))) {
    meta <- .sc_sample_metadata_row(sample_metadata, ids[[i]])
    read_args <- list(
      path = found$matrix_path[[i]], sample_id = ids[[i]], format = found$format[[i]],
      version = version, data_type = data_type,
      assay = assay, project = ids[[i]], metadata = meta,
      delimiter = delimiter, return_report = TRUE,
      na_to_zero = na_to_zero
    )
    if (!orientation_was_missing) read_args$orientation <- orientation
    ans <- do.call(read_sc_data, read_args)
    objects[[i]] <- .sc_set_sample_identity(ans$object, ids[[i]], set_orig_ident = TRUE)
    reports[[i]] <- ans$report
  }
  if (!merge) {
    if (return_report) return(list(objects = objects, report = reports))
    return(objects)
  }
  merged <- if (length(objects) == 1L) objects[[1L]] else {
    base::merge(
      x = objects[[1L]], y = objects[-1L],
      add.cell.ids = names(objects),
      project = basename(normalizePath(path, winslash = "/", mustWork = TRUE)),
      merge.data = identical(data_type, "data")
    )
  }
  if (return_report) return(list(object = merged, report = reports))
  merged
}
