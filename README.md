<div align="center">
  <h1>📜 <code>immunread</code> – <strong>Read immune repertoire data in R</strong></h1>
</div>

---

<div align="center">
  <a href="https://github.com/immunomind">
    <img alt="Ecosystem: ImmunoMind"
         src="https://img.shields.io/badge/ecosystem-ImmunoMind-orange?style=flat-square">
  </a>
  <a href="https://github.com/immunomind/immunread/actions/workflows/altdoc.yaml">
    <img alt="Build documentation"
         src="https://github.com/immunomind/immunread/actions/workflows/altdoc.yaml/badge.svg">
  </a>
  <a href="https://github.com/immunomind/immunread/issues">
    <img alt="GitHub Issues"
         src="https://img.shields.io/github/issues/immunomind/immunread?style=flat-square">
  </a>
</div>

<p align="center">
  <a href="https://immunomind.github.io/immunread/">Documentation</a>
  |
  <a href="https://immunomind.github.io/docs/">Ecosystem</a>
</p>

---

`immunread` simplifies reading immune repertoire data from the most common
bulk and single-cell formats into `immundata::ImmunData` objects.

## Supported parsers

| Tool / standard | Input format | Parser |
|---|---|---|
| AIRR-C | Bulk Rearrangement table (`.tsv`) | [`read_airr_bulk()`](https://immunomind.github.io/immunread/man/read_airr_bulk.html) |
| | Single-cell Rearrangement table (`.tsv`) | [`read_airr_sc()`](https://immunomind.github.io/immunread/man/read_airr_sc.html) |
| 10x Genomics | V(D)J contig annotations (`filtered_contig_annotations.csv`, `all_contig_annotations.csv`) | [`read_10x()`](https://immunomind.github.io/immunread/man/read_10x.html) |
| Adaptive immunoSEQ | Analyzer export | [`read_adaptive()`](https://immunomind.github.io/immunread/man/read_adaptive.html) |
| | Analyzer v2 export | [`read_adaptive()`](https://immunomind.github.io/immunread/man/read_adaptive.html) |
| | Connector export | [`read_adaptive()`](https://immunomind.github.io/immunread/man/read_adaptive.html) |
| BD Rhapsody | VDJ Contigs AIRR table (`*_VDJ_Dominant_Contigs_AIRR.tsv`, `*_VDJ_Unfiltered_Contigs_AIRR.tsv`) | [`read_bd()`](https://immunomind.github.io/immunread/man/read_bd.html) |
| Dandelion | AIRR contig table (`Dandelion.data`) | [`read_dandelion()`](https://immunomind.github.io/immunread/man/read_dandelion.html) |
| Immcantation | AIRR-C bulk table (`.tsv`) | [`read_immcantation_bulk()`](https://immunomind.github.io/immunread/man/read_immcantation_bulk.html) |
| | AIRR-C single-cell table (`.tsv`) | [`read_immcantation_sc()`](https://immunomind.github.io/immunread/man/read_immcantation_sc.html) |
| MiXCR | AIRR-C bulk export (`.tsv`) | [`read_airr_bulk()`](https://immunomind.github.io/immunread/man/read_airr_bulk.html) |
| | AIRR-C single-cell export (`.tsv`) | [`read_airr_sc()`](https://immunomind.github.io/immunread/man/read_airr_sc.html) |
| | Single-cell clone export | [`read_mixcr()`](https://immunomind.github.io/immunread/man/read_mixcr.html) |
| Parse Biosciences | Wide per-cell TCR table | [`read_parsebio()`](https://immunomind.github.io/immunread/man/read_parsebio.html) |
| | AIRR-C export (`tcr_annotation_airr.tsv`) | [`read_airr_sc()`](https://immunomind.github.io/immunread/man/read_airr_sc.html) |
| TRUST4 | Cell table with embedded chain records | [`read_trust4()`](https://immunomind.github.io/immunread/man/read_trust4.html) |
| WAT3R | Per-cell TCR table (`barcode_results.csv`) | [`read_wat3r()`](https://immunomind.github.io/immunread/man/read_wat3r.html) |

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
