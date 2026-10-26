#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(scAutoLoad))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
package_root <- if (length(script_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/")))
} else {
  normalizePath(getwd(), winslash = "/")
}
case_root <- file.path(package_root, "demo_cases")

cases <- list(
  mtx_standard = list(
    path = "标准10X读入------.GSE213216", format = "10x_mtx", batch = TRUE,
    description = "10x MTX triplets in per-sample folders"
  ),
  h5_single = list(
    path = file.path("单个h5文件", "GSM5861295.h5"), format = "10x_h5", batch = FALSE,
    description = "One 10x HDF5 file"
  ),
  rds_saved_object = list(
    path = file.path("单个h5文件", "af.rds"), format = "rds", batch = FALSE,
    description = "A saved Seurat object in RDS format"
  ),
  txt_multi = list(
    path = "多个矩阵读入", format = "txt", batch = FALSE, sequential = TRUE,
    description = "Several tab-delimited gene-by-cell matrices"
  ),
  csv_gz_transposed = list(
    path = "多个gz倒置xls------GSE118828", format = "csv", batch = TRUE,
    description = "Gzipped CSV with cells in rows and genes in columns"
  ),
  csv_gz_multi = list(
    path = "多个gz的csv文件", format = "csv", batch = TRUE,
    description = "Several gzipped CSV count matrices"
  ),
  tsv_gz_multi = list(
    path = "多个gz的tsv文件", format = "tsv", batch = TRUE,
    description = "Several gzipped TSV count matrices"
  ),
  h5_multi = list(
    path = "多个H5读入", format = "10x_h5", batch = TRUE,
    description = "Several 10x HDF5 files"
  ),
  geo_text = list(
    path = file.path("GEO有分号读入（单个矩阵）--课程示例数据",
                     "GSE146115_HCC1-2-5-9_count_with_ERCC.txt"),
    format = "txt", batch = FALSE,
    description = "GEO text matrix; delimiter is detected from file contents"
  ),
  csv_gz_single = list(
    path = file.path("gz的单个csv文件", "GSE115469_Data.csv.gz"),
    format = "csv", batch = FALSE,
    description = "One gzipped CSV count matrix"
  )
)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || identical(args[[1L]], "list")) {
  cat("Available cases:\n")
  for (nm in names(cases)) cat("  ", nm, " - ", cases[[nm]]$description, "\n", sep = "")
  cat("\nRun one case at a time, for example:\n  Rscript examples/run_real_cases.R mtx_standard\n")
  quit(status = 0L)
}

case_name <- args[[1L]]
if (!case_name %in% names(cases)) {
  stop("Unknown case: ", case_name, ". Run this script with 'list' to see choices.")
}
if (!dir.exists(case_root)) {
  stop("Could not find demo_cases at ", case_root,
       ". Run this script from the scAutoLoad source checkout.")
}

case <- cases[[case_name]]
input_path <- file.path(case_root, case$path)
cat("Case:", case_name, "-", case$description, "\n")
print(scan_sc_data(input_path))

if (case$batch) {
  objects <- read_sc_batch(
    input_path,
    formats = case$format,
    data_type = "counts",
    version = "auto",
    merge = FALSE
  )
  dims <- do.call(rbind, lapply(objects, function(x) dim(x)))
  colnames(dims) <- c("features", "cells")
  print(dims)
  failures <- vapply(objects, function(x) {
    sum(validate_sc_data(x)$status == "fail")
  }, integer(1))
  print(data.frame(sample_id = names(objects), validation_failures = failures,
                   row.names = NULL))
} else if (isTRUE(case$sequential)) {
  candidates <- scan_sc_data(input_path)
  candidates <- candidates[candidates$format %in% case$format, , drop = FALSE]
  if (!nrow(candidates)) stop("No matching input files were found for this case.")
  for (i in seq_len(nrow(candidates))) {
    cat("\nSample:", candidates$sample_id[[i]], "\n")
    object <- read_sc_data(
      candidates$matrix_path[[i]],
      data_type = "counts",
      version = "auto"
    )
    print(dim(object))
    print(validate_sc_data(object))
    rm(object)
    invisible(gc())
  }
} else {
  object <- read_sc_data(
    input_path,
    data_type = "counts",
    version = "auto"
  )
  print(dim(object))
  print(validate_sc_data(object))
}

# To force a known orientation for a single matrix, pass orientation = TRUE
# (features in rows) or orientation = FALSE (features in columns) to read_sc_data().

