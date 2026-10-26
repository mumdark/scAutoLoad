\name{read_sc_batch}
\alias{read_sc_batch}
\title{Read All Detected Samples in a Directory}
\usage{
read_sc_batch(path, recursive = TRUE, sample_metadata = NULL, merge = FALSE,
  version = "auto", orientation = TRUE, data_type = "counts",
  assay = "RNA", delimiter = NULL, return_report = FALSE, formats = NULL,
  na_to_zero = TRUE)
}
\arguments{
\item{path}{Directory to scan.}
\item{recursive}{Whether to scan nested directories.}
\item{sample_metadata}{Optional sample-level data frame keyed by \code{sample_id}.}
\item{merge}{Combine samples into one Seurat object. Sample names are prefixed to cell barcodes and stored as \code{orig.ident}. For \code{data_type = "data"}, normalized values are retained during the merge.}
\item{version}{Seurat assay layout, passed to \code{read_sc_data()}.}
\item{orientation}{Logical orientation: \code{TRUE} means features in rows and cells in columns; \code{FALSE} means cells in rows and features in columns. If omitted, each sample is auto-detected and ambiguous matrices use \code{TRUE}.}
\item{data_type}{\code{"counts"} by default; use \code{"data"} for normalized expression.}
\item{assay}{Assay name.}
\item{delimiter}{Optional delimiter override for text files.}
\item{return_report}{If \code{TRUE}, return imported object(s) with a combined report.}
\item{formats}{Optional character vector selecting formats from the scan, such as \code{"10x_mtx"}, \code{"csv"}, \code{"tsv"}, \code{"txt"}, \code{"10x_h5"}, or \code{"rds"}.}
\item{na_to_zero}{If \code{TRUE} (default), replace missing expression values with zero and
print how many values were replaced for each sample. Set \code{FALSE} to retain them and
stop before creating the Seurat object.}
}
\value{
A named list of Seurat objects, a merged Seurat object, or a list containing
object(s) and the import report.
}
\description{
Scans a directory and reads each detected sample. Incomplete 10x MTX datasets
and duplicate sample IDs are reported as errors rather than silently skipped.
}

