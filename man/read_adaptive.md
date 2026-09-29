

# Read Adaptive immunoSEQ Analyzer and Connector tables

## Description

Reads Adaptive immunoSEQ Analyzer or Connector tables and maps their
resolved gene, CDR3, sequence, and abundance fields to AIRR-style
columns.

Resolved gene calls are used directly. <code>templates</code> supplies
bulk abundance for Analyzer and Connector; Analyzer v2 uses <code>count
(templates/reads)</code>. <code>seq_reads</code> remains a distinct
source count when present.

## Usage

<pre><code class='language-R'>read_adaptive(
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

path <- system.file("extdata", "adaptive", "analyzer", "contigs.tsv.gz",
  package = "immunread", mustWork = TRUE)
idata <- read_adaptive(path, schema = c("v_call", "junction_aa"),
  output_folder = tempfile("immunread-"))
```
