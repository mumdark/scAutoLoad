.sc_validation_row <- function(check, status, details) {
  data.frame(check = check, status = status, details = details, stringsAsFactors = FALSE)
}

#' Validate an imported single-cell object
#'
#' @param object Seurat object to inspect.
#' @param deep Whether to inspect expression values for negative, non-finite,
#'   and non-integer values. Deep checks may take time for large matrices.
#' @return A data frame with one row per validation check.
#' @export
validate_sc_data <- function(object, deep = FALSE) {
  if (!inherits(object, "Seurat")) .sc_abort("`object` must be a Seurat object.")
  .sc_require_seuratobject()
  assay_name <- SeuratObject::DefaultAssay(object)
  assay_obj <- object[[assay_name]]
  rn <- rownames(assay_obj)
  cn <- colnames(assay_obj)
  meta <- object@meta.data

  layer_names <- if (exists("Layers", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character())
  } else {
    count_mat <- tryCatch(SeuratObject::GetAssayData(object, assay = assay_name, slot = "counts"),
                          error = function(e) NULL)
    data_mat <- tryCatch(SeuratObject::GetAssayData(object, assay = assay_name, slot = "data"),
                         error = function(e) NULL)
    out <- character()
    if (!is.null(count_mat) && length(dim(count_mat)) == 2L && all(dim(count_mat) > 0L)) out <- c(out, "counts")
    if (!is.null(data_mat) && length(dim(data_mat)) == 2L && all(dim(data_mat) > 0L)) out <- c(out, "data")
    out
  }
  count_layers <- layer_names[grepl("^counts($|\\.)", layer_names)]
  data_layers <- layer_names[grepl("^data($|\\.)", layer_names)]
  selected_layer <- if (length(count_layers)) count_layers[[1L]] else if (length(data_layers)) data_layers[[1L]] else NA_character_
  supported_assay <- inherits(assay_obj, "Assay5") || inherits(assay_obj, "Assay")

  checks <- list(
    .sc_validation_row("seurat_class", "pass", paste(class(object), collapse = "/")),
    .sc_validation_row("non_empty_dimensions", if (length(rn) > 0L && length(cn) > 0L) "pass" else "fail",
                       paste(length(rn), "features x", length(cn), "cells")),
    .sc_validation_row("unique_feature_names", if (!anyDuplicated(rn)) "pass" else "fail",
                       if (!anyDuplicated(rn)) "Feature names are unique." else "Feature names contain duplicates."),
    .sc_validation_row("unique_cell_names", if (!anyDuplicated(cn)) "pass" else "fail",
                       if (!anyDuplicated(cn)) "Cell names are unique." else "Cell names contain duplicates."),
    .sc_validation_row("metadata_alignment", if (identical(rownames(meta), colnames(object))) "pass" else "fail",
                       "Metadata row names should match object cell names in order."),
    .sc_validation_row("expression_layer", if (!is.na(selected_layer)) "pass" else "fail",
                       if (is.na(selected_layer)) "No counts or data layer found." else paste("Using layer", selected_layer)),
    .sc_validation_row("counts_layer", if (length(count_layers)) "pass" else "warning",
                       if (length(count_layers)) paste(count_layers, collapse = ", ") else "No raw counts layer; object may contain normalized data only."),
    .sc_validation_row("assay_layout", if (supported_assay) "pass" else "fail",
                       paste(class(assay_obj), collapse = "/"))
  )

  if (deep && !is.na(selected_layer)) {
    mat <- tryCatch({
      if (inherits(assay_obj, "Assay5")) {
        SeuratObject::LayerData(object, assay = assay_name, layer = selected_layer)
      } else if (length(count_layers)) {
        SeuratObject::GetAssayData(object, assay = assay_name, slot = "counts")
      } else {
        SeuratObject::GetAssayData(object, assay = assay_name, slot = "data")
      }
    }, error = function(e) e)
    if (inherits(mat, "error")) {
      checks[[length(checks) + 1L]] <- .sc_validation_row("expression_values", "warning", conditionMessage(mat))
    } else {
      flags <- .sc_value_flags(mat)
      finite <- flags[["finite"]]
      nonnegative <- flags[["nonnegative"]]
      integerish <- flags[["integerish"]]
      checks[[length(checks) + 1L]] <- .sc_validation_row(
        "finite_values", if (finite) "pass" else "fail",
        if (finite) "All inspected values are finite." else "NA, NaN, or infinite values were found."
      )
      checks[[length(checks) + 1L]] <- .sc_validation_row(
        "nonnegative_values", if (nonnegative) "pass" else "warning",
        if (nonnegative) "All inspected values are nonnegative." else "Negative values are not valid raw counts."
      )
      checks[[length(checks) + 1L]] <- .sc_validation_row(
        "integer_counts", if (integerish) "pass" else "warning",
        if (integerish) "Values are integer-like." else "Values include decimals; confirm the selected layer represents raw counts."
      )
    }
  }

  out <- do.call(rbind, checks)
  rownames(out) <- NULL
  class(out) <- c("sc_validation", class(out))
  out
}

