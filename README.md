<div align="center">
  <h1>📜 <code>immunread</code> — Read immune repertoire data in R</h1>
</div>

<p align="center">
  <a href="https://immunomind.github.io/immunread/">API reference</a>
  · <a href="https://github.com/immunomind/immunread/issues">Issues</a>
  · <a href="https://github.com/immunomind">ImmunoMind ecosystem</a>
</p>

`immunread` imports
[AIRR](https://docs.airr-community.org/en/stable/datarep/rearrangements.html),
[10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/cr-5p-outputs-annotations-vdj),
[BD Rhapsody](https://bd-rhapsody-bioinfo-docs.genomics.bd.com/outputs/outputs_vdj_contigs.html),
[MiXCR](https://mixcr.com/mixcr/reference/overview-analysis-overview/),
[Immcantation](https://immcantation.readthedocs.io/en/stable/),
[Dandelion](https://sc-dandelion.readthedocs.io/en/latest/index.html),
[Parse Biosciences](https://www.parsebiosciences.com/products/evercode-tcr/),
[TRUST4](https://github.com/liulab-dfci/TRUST4),
[WAT3R](https://github.com/mainciburu/WAT3R), and
[Adaptive immunoSEQ](https://www.adaptivebiotech.com/) data into
[`immundata::ImmunData`](https://github.com/immunomind/immundata) objects.

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
