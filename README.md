<div align="center">
  <h1>📜 <code>immunread</code> — Read immune repertoire data in R</h1>
</div>

<p align="center">
  <a href="https://immunomind.github.io/immunread/">API reference</a>
  · <a href="https://github.com/immunomind/immunread/issues">Issues</a>
  · <a href="https://github.com/immunomind">ImmunoMind ecosystem</a>
</p>

`immunread` simplifies reading immune repertoire data from the most common
bulk and single-cell formats into `immundata::ImmunData` objects.

## Supported parsers

- AIRR Rearrangement tables: [bulk](https://immunomind.github.io/immunread/man/read_airr_bulk.html) and [single-cell](https://immunomind.github.io/immunread/man/read_airr_sc.html)
- [10x Genomics](https://immunomind.github.io/immunread/man/read_10x.html)
- [Adaptive immunoSEQ](https://immunomind.github.io/immunread/man/read_adaptive.html)
- [BD Rhapsody](https://immunomind.github.io/immunread/man/read_bd.html)
- [Dandelion](https://immunomind.github.io/immunread/man/read_dandelion.html)
- Immcantation: [bulk](https://immunomind.github.io/immunread/man/read_immcantation_bulk.html) and [single-cell](https://immunomind.github.io/immunread/man/read_immcantation_sc.html)
- [MiXCR](https://immunomind.github.io/immunread/man/read_mixcr.html)
- [Parse Biosciences](https://immunomind.github.io/immunread/man/read_parsebio.html)
- [TRUST4](https://immunomind.github.io/immunread/man/read_trust4.html)
- [WAT3R](https://immunomind.github.io/immunread/man/read_wat3r.html)

## Installation

Install the development packages from GitHub with
[`pak`](https://pak.r-lib.org/) in `R`:

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

The package is under development, so its API may change before version 1.0.
