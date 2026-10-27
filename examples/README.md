# scAutoLoad: Real-Data Examples

This directory contains examples for the real single-cell input cases. The
datasets occupy about 1.8 GiB and are not included in the GitHub repository or
R package builds. Run these examples from the package source directory after
copying the case folders into `demo_cases/`.

## Installation and parameters

```r
install.packages("remotes")
remotes::install_github("mumdark/scAutoLoad")
```

Reading 10x HDF5 files requires `Seurat`; creating Seurat objects requires
`SeuratObject`. `orientation` is logical. If omitted, scAutoLoad detects the
matrix direction from identifiers and falls back to `TRUE` (features in rows,
cells in columns) when ambiguous. Pass `TRUE` or `FALSE` to force a direction.
Matrix files are treated as counts by default.

`na_to_zero = TRUE` is the default. Missing expression values, including blank
fields, are replaced by zero and the number replaced is printed. Set
`na_to_zero = FALSE` to preserve missing values and stop with an error. This
option applies to expression matrices and does not modify saved Seurat RDS
objects.

## R calls for every case

Run only the calls you need. Multi-sample calls return one merged Seurat object
by default. Sample IDs are inferred from the input file names or sample-folder
names and are added to the `sample_id` and `orig.ident` metadata columns. Pass
`merge = FALSE` to keep a named list of per-sample objects. For the large TXT
case, use the sequential runner below to limit memory use.

```r
library(scAutoLoad)
root <- "demo_cases"

# 1. Standard 10x MTX: AD1, AD2, CO1, and CO2 samples
mtx <- read_sc_batch(
  file.path(root, "标准10X读入------.GSE213216"),
  formats = "10x_mtx", version = "auto"
)

# 2. One 10x HDF5 file
h5_one <- read_sc_data(
  file.path(root, "单个h5文件", "GSM5861295.h5"), version = "auto"
)

# 3. A saved Seurat object in RDS format
rds_object <- read_sc_data(
  file.path(root, "单个h5文件", "af.rds"), version = "auto"
)

# 4. Multiple plain TXT matrices
# Use the runner below to read D1, D2, and D3 sequentially and release each.
# txt_objects <- read_sc_batch(
#   file.path(root, "多个矩阵读入"), formats = "txt", merge = FALSE
# )

# 5. Transposed gzip CSV files: cells in rows, features in columns
csv_transposed <- read_sc_batch(
  file.path(root, "多个gz倒置xls------GSE118828"),
  formats = "csv", version = "auto"
)

# 6. Multiple .csv.gz matrices
csv_gz <- read_sc_batch(
  file.path(root, "多个gz的csv文件"),
  formats = "csv", version = "auto"
)

# 7. Multiple .tsv.gz matrices
tsv_gz <- read_sc_batch(
  file.path(root, "多个gz的tsv文件"),
  formats = "tsv", version = "auto"
)

# 8. Multiple 10x HDF5 files
h5_many <- read_sc_batch(
  file.path(root, "多个H5读入"),
  formats = "10x_h5", version = "auto"
)

# 9. GEO TXT matrix; the delimiter is TAB and blank values are filled with 0
geo_txt <- read_sc_data(
  file.path(root, "GEO有分号读入（单个矩阵）--课程示例数据",
            "GSE146115_HCC1-2-5-9_count_with_ERCC.txt"),
  delimiter = "\t", na_to_zero = TRUE, version = "auto"
)

# 10. One .csv.gz matrix
csv_gz_one <- read_sc_data(
  file.path(root, "gz的单个csv文件", "GSE115469_Data.csv.gz"),
  version = "auto"
)

# Inspect dimensions and validity, for example:
dim(csv_transposed)
validate_sc_data(csv_transposed)

# Force cells in rows and features in columns for a known input layout:
# read_sc_data("sample.csv.gz", orientation = FALSE)

# Keep missing values as an error instead of converting them to zero:
# read_sc_data("sample.csv.gz", na_to_zero = FALSE)
```

## Run each integration case

From a terminal in the package source directory, run one case at a time:

```sh
Rscript examples/run_real_cases.R list
Rscript examples/run_real_cases.R mtx_standard
Rscript examples/run_real_cases.R h5_single
Rscript examples/run_real_cases.R rds_saved_object
Rscript examples/run_real_cases.R txt_multi
Rscript examples/run_real_cases.R csv_gz_transposed
Rscript examples/run_real_cases.R csv_gz_multi
Rscript examples/run_real_cases.R tsv_gz_multi
Rscript examples/run_real_cases.R h5_multi
Rscript examples/run_real_cases.R geo_text
Rscript examples/run_real_cases.R csv_gz_single
```

The runner prints the scan results, imported dimensions, and validation checks.
It reads the three large TXT matrices one at a time and releases each object
before continuing.
