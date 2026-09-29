

# Read TRUST4 cell tables

## Description

Reads TRUST4 cell tables, expands embedded chain records, and maps their
receptor and count fields to AIRR-style columns.

Each supplied primary or secondary chain becomes a row.
<code>immundata</code> ranks chains when constructing receptors. The
embedded <code>read_cnt</code> becomes <code>consensus_count</code>,
following TRUST4’s AIRR export.

## Usage

<pre><code class='language-R'>read_trust4(
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

path <- system.file("extdata", "trust4", "1.1.7", "contigs.csv.gz",
  package = "immunread", mustWork = TRUE)
idata <- read_trust4(path, schema = c("v_call", "junction_aa"),
  output_folder = tempfile("immunread-"))
```
