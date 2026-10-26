\name{read_sc_data}
\alias{read_sc_data}
\title{Read a Single-Cell Dataset into a Seurat Object}
\usage{
read_sc_data(path, sample_id = NULL, format = "auto", version = "auto",
  orientation = TRUE, data_type = "counts", assay = "RNA",
  project = NULL, metadata = NULL, delimiter = NULL, recursive = TRUE,
  return_report = FALSE, na_to_zero = TRUE)
}
\arguments{
\item{path}{A matrix file, RDS/HDF5 file, or directory containing one dataset.}
\item{sample_id}{Optional sample identifier; defaults to a name inferred from the path.}
\item{format}{Input format, or \code{"auto"} to detect it.}
\item{version}{Seurat assay layout: \code{"auto"}, \code{"v4"} (legacy Assay), or \code{"v5"} (Assay5).}
\item{orientation}{Logical orientation: \code{TRUE} means features in rows and cells in columns; \code{FALSE} means cells in rows and features in columns. If omitted, orientation is detected from identifiers and falls back to \code{TRUE} when ambiguous. Explicitly passing \code{TRUE} or \code{FALSE} forces it.}
\item{data_type}{\code{"counts"} by default; use \code{"data"} for normalized expression.}
\item{na_to_zero}{If \code{TRUE} (default), replace missing expression values with zero and
print how many values were replaced. Set \code{FALSE} to keep them and stop with an error
before creating the Seurat object.}
\item{assay}{Assay name for the expression matrix.}
\item{project}{Project name stored in the Seurat object.}
\item{metadata}{Optional cell-level data frame or delimited file. Cell names must be row names or provided in a \code{cell} column; a single row is treated as sample-level metadata.}
\item{delimiter}{Optional delimiter override for text matrices: comma, tab, or semicolon.}
\item{recursive}{Scan nested directories when \code{path} is a directory.}
\item{return_report}{If \code{TRUE}, return \code{list(object, report)}.}
}
\value{
A Seurat object, or an object and import report when requested.
}
\description{
Reads supported single-cell formats, detects or applies matrix orientation,
and creates an object using the requested Seurat assay layout.
}

