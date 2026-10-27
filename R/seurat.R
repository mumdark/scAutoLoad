.sc_require_seuratobject <- function() {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    .sc_abort("Creating Seurat objects requires 'SeuratObject'. Install it with install.packages('SeuratObject').")
  }
}

.sc_resolve_version <- function(version = "auto") {
  .sc_require_seuratobject()
  version <- match.arg(version, c("auto", "v4", "v5"))
  installed <- utils::packageVersion("SeuratObject")
  if (version == "auto") {
    major <- as.integer(strsplit(as.character(installed), "\\.")[[1L]][[1L]])
    if (!is.na(major) && major >= 5L) return("v5")
    return("v4")
  }
  if (version == "v5" && !exists("CreateAssay5Object", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    .sc_abort("`version = 'v5'` requires SeuratObject 5 or later; installed version is ", installed, ".")
  }
  version
}

.sc_create_seurat <- function(mat, data_type = "counts", assay = "RNA", project = "scAutoLoad",
                             metadata = NULL, version = "auto") {
  target <- .sc_resolve_version(version)
  mat <- .sc_to_sparse(mat)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    .sc_abort("The expression matrix must have feature names and cell names before creating a Seurat object.")
  }
  args <- if (data_type == "data") list(data = mat) else list(counts = mat)
  if (target == "v5") {
    assay_obj <- do.call(SeuratObject::CreateAssay5Object, args)
  } else {
    assay_obj <- do.call(SeuratObject::CreateAssayObject, args)
  }
  object <- SeuratObject::CreateSeuratObject(
    counts = assay_obj, assay = assay, project = project, meta.data = metadata
  )
  list(object = object, assay_version = target,
       seuratobject_version = as.character(utils::packageVersion("SeuratObject")))
}

.sc_set_sample_identity <- function(object, sample_id, set_orig_ident = FALSE) {
  metadata <- object@meta.data
  metadata$sample_id <- as.character(sample_id)[1L]
  if (isTRUE(set_orig_ident)) metadata$orig.ident <- as.character(sample_id)[1L]
  object@meta.data <- metadata
  object
}

.sc_read_metadata <- function(metadata, cell_names, sample_id) {
  if (is.null(metadata)) return(NULL)
  if (is.character(metadata) && length(metadata) == 1L) {
    if (!file.exists(metadata)) .sc_abort("Metadata file does not exist: ", metadata)
    first <- .sc_read_lines(metadata, n = 1L)
    sep <- .sc_detect_separator(if (length(first)) first[[1L]] else "")
    con <- .sc_open_text(metadata)
    on.exit(close(con), add = TRUE)
    dat <- utils::read.table(con, header = TRUE, sep = sep, check.names = FALSE,
                             stringsAsFactors = FALSE, comment.char = "")
    metadata <- dat
  }
  if (!is.data.frame(metadata)) .sc_abort("`metadata` must be a data frame or a delimited-file path.")
  implicit_rownames <- identical(rownames(metadata), as.character(seq_len(nrow(metadata))))
  if ("cell" %in% tolower(names(metadata))) {
    cell_col <- which(tolower(names(metadata)) == "cell")[[1L]]
    rownames(metadata) <- as.character(metadata[[cell_col]])
    if (!all(cell_names %in% rownames(metadata))) {
      .sc_abort("Metadata 'cell' values do not match all expression-matrix cell names.")
    }
    metadata <- metadata[cell_names, , drop = FALSE]
    metadata[[cell_col]] <- NULL
  } else if (!implicit_rownames && all(cell_names %in% rownames(metadata))) {
    metadata <- metadata[cell_names, , drop = FALSE]
  } else if (nrow(metadata) == 1L && implicit_rownames) {
    metadata <- metadata[rep(1L, length(cell_names)), , drop = FALSE]
    rownames(metadata) <- cell_names
  } else {
    .sc_abort("Metadata must have cell names as row names, a 'cell' column, or exactly one sample-level row.")
  }
  metadata$orig.ident <- sample_id
  metadata
}
