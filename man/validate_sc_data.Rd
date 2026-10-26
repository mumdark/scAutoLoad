\name{validate_sc_data}
\alias{validate_sc_data}
\title{Validate an Imported Single-Cell Object}
\usage{
validate_sc_data(object, deep = FALSE)
}
\arguments{
\item{object}{Seurat object to inspect.}
\item{deep}{Whether to inspect expression values for negative, non-finite, and non-integer values. Deep checks may take time for large matrices.}
}
\value{
A data frame with one row per validation check.
}
\description{
Checks Seurat class, dimensions, unique names, metadata alignment, and the
presence of counts or data layers.
}

