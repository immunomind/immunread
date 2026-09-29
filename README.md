<div align="center">
  <h1>📥 <code>immunread</code> — Read immune repertoire data in R</h1>
</div>

<p align="center">
  <a href="https://immunomind.github.io/immunread/">API reference</a>
  · <a href="https://github.com/immunomind/immunread/issues">Issues</a>
  · <a href="https://github.com/immunomind">ImmunoMind ecosystem</a>
</p>

`immunread` imports AIRR, 10x Genomics, BD Rhapsody, MiXCR, Immcantation,
Dandelion, Parse Biosciences, TRUST4, WAT3R, and Adaptive immunoSEQ data into
[`immundata::ImmunData`](https://github.com/immunomind/immundata) objects.
Readers use duckplyr for large datasets and accept manifests through
`immundata::read_repertoires()`.

## Installation

Install the development packages from GitHub:

```r
install.packages("pak")
pak::pak(c("immunomind/immundata@dev", "immunomind/immunread"))
```

## Example

```r
library(immunread)

path <- system.file(
  "extdata", "airr", "2.0-bulk", "contigs.tsv.gz",
  package = "immunread"
)
idata <- read_airr_bulk(
  path,
  schema = c("v_call", "junction_aa"),
  output_folder = tempfile("immunread-")
)
```

The package is under development; its API may change before version 1.0.
