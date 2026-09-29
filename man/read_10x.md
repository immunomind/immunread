

# Read 10x Genomics V(D)J contigs

## Description

Reads 10x Genomics V(D)J contig annotations and maps their barcode,
locus, gene, sequence, and count fields to AIRR-style columns.

Supplied secondary chains are retained; <code>immundata</code> ranks
chains when constructing receptors. <code>umis</code> and
<code>reads</code> remain distinct counts.

## Usage

<pre><code class='language-R'>read_10x(
  path,
  schema,
  manifest = NULL,
  manifest_file_col = "file",
  column_profile = "standard",
  repertoire_schema = "&lt;auto&gt;",
  ...
)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="path">path</code>
</td>
<td>
Repertoire file paths, glob, or <code>“\<manifest\>”</code>.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="schema">schema</code>
</td>
<td>
Receptor features or an <code>immundata::make_receptor_schema()</code>
result.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="manifest">manifest</code>
</td>
<td>
Optional native immundata manifest.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="manifest_file_col">manifest_file_col</code>
</td>
<td>
Manifest column containing file paths when <code>path =
“\<manifest\>”</code>.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="column_profile">column_profile</code>
</td>
<td>
Source columns to keep: <code>“minimal”</code>, <code>“standard”</code>,
or <code>“full”</code>.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="repertoire_schema">repertoire_schema</code>
</td>
<td>
Repertoire grouping passed to
<code>immundata::read_repertoires()</code>.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="...">…</code>
</td>
<td>
Other arguments passed to <code>immundata::read_repertoires()</code>.
</td>
</tr>
</table>

## Examples

``` r
library("immunread")

path <- system.file("extdata", "10x", "9.0.1", "contigs.csv.gz",
  package = "immunread", mustWork = TRUE)
idata <- read_10x(path, schema = c("v_call", "junction_aa"),
  output_folder = tempfile("immunread-"))
```
