# scAutoLoad

`scAutoLoad` detects common single-cell expression matrix formats, reads one or
many samples, reports matrix orientation, and creates Seurat objects compatible
with the installed `SeuratObject` version.

## Supported formats

- 10x Matrix Market triplets: `matrix.mtx[.gz]`, `features/genes.tsv[.gz]`,
  and `barcodes.tsv[.gz]`
- Dense CSV, TSV, and TXT matrices, including gzip-compressed files and
  nonstandard names such as `.gz.csv`, detected using the gzip signature
- 10x Cell Ranger HDF5 files (`.h5`; requires `Seurat`)
- RDS files containing a Seurat object, a `SingleCellExperiment`, or a numeric
  matrix
- Single-sample and multi-sample directory layouts
- Seurat v4 `Assay` and Seurat v5 `Assay5` output, selected from the installed
  `SeuratObject` version

AnnData/H5AD and Loom are planned for a later release.

## Installation

```r
install.packages(c("Matrix", "SeuratObject"))

# Required only for reading 10x HDF5 files.
install.packages("Seurat")

# Optional: faster parsing for delimited text matrices.
install.packages("data.table")

install.packages("remotes")
remotes::install_github("mumdark/scAutoLoad")
```

To install a local checkout, use `remotes::install_local("path/to/scAutoLoad")`.

## Scan and read

```r
library(scAutoLoad)

scan_sc_data("./GEO_data")

obj <- read_sc_data(
  "./GEO_data/sample1",
  data_type = "counts", # default for matrix files
  version = "auto"
)

obj_merged <- read_sc_batch("./GEO_data") # merge = TRUE by default
table(obj_merged$sample_id)

# Keep samples as separate objects when needed:
objects <- read_sc_batch("./GEO_data", merge = FALSE)
validate_sc_data(obj)
```

`orientation` is logical: `TRUE` means features in rows and cells in columns;
`FALSE` means cells in rows and features in columns. If omitted,
`read_sc_data()` detects direction from identifiers, prints its decision, and
uses `TRUE` when the direction is ambiguous. Explicitly pass `TRUE` or `FALSE`
to force an orientation.

Missing expression values are replaced with zero by default
(`na_to_zero = TRUE`), and the number replaced is printed. Set
`na_to_zero = FALSE` to retain missing values and stop with an error. Counts are
the default input type. Non-integer values produce a warning so that
normalized expression matrices are not mistaken for raw counts.

Duplicate feature names are made unique without summing rows; the original
names are retained in the import report. `read_sc_batch()` merges all detected
samples into one Seurat object by default. It infers each sample ID from its
file name or containing folder, prefixes cell names with that ID, and records
it in both the `sample_id` and `orig.ident` metadata columns. Set `merge = FALSE`
to return one object per sample.

Delimited text matrices are read into memory before conversion to a sparse
matrix. Large dense CSV/TSV/TXT files can require several times their file size
in RAM; 10x MTX is preferred for large sparse count matrices.

## Real-data examples

The `examples/` directory contains R calls for each tested input layout and a
runner that processes one case at a time. The original case data is about
1.8 GiB and is not included in the GitHub repository or package builds. To run
the integration examples, copy the supplied case folders into `demo_cases/`
under the source checkout, then follow [`examples/README.md`](examples/README.md).
