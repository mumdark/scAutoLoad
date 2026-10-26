\name{scan_sc_data}
\alias{scan_sc_data}
\title{Scan a Directory for Single-Cell Datasets}
\usage{
scan_sc_data(path, recursive = TRUE,
  orientation = TRUE)
}
\arguments{
\item{path}{A file or directory to inspect.}
\item{recursive}{Whether to scan nested directories.}
\item{orientation}{Logical fallback orientation: \code{TRUE} means features in rows and cells in columns; \code{FALSE} means cells in rows and features in columns.}
}
\value{
A data frame with one row per candidate dataset and class \code{sc_scan}.
}
\description{
Inspects file names, gzip signatures, table headers, and 10x sidecar files.
Scanning does not read complete expression matrices.
}

