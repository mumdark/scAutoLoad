\name{print.sc_scan}
\alias{print.sc_scan}
\title{Print a Single-Cell Scan Report}
\usage{
\method{print}{sc_scan}(x, ...)
}
\arguments{
\item{x}{An object returned by \code{scan_sc_data()}.}
\item{...}{Additional arguments passed to \code{print.data.frame()}.}
}
\value{
The input scan report, invisibly.
}
\description{
Displays detected sample IDs, formats, orientation, status, and notes.
}

