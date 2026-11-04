

# Read MiXCR single-cell clones

[**Source code**](https://github.com/immunomind/immunread/tree/main/R/read_mixcr.R#L41)

## Description

Reads MiXCR
<a href="https://mixcr.com/mixcr/reference/mixcr-export/#single-cell-libraries">single-cell
clone exports</a> and maps their cell, chain, gene, sequence, and count
fields to AIRR-style columns.

The input must have <code>cellGroup</code>, <code>topChains</code>, and
<code>uniqueMoleculeCount</code> columns, as in the bundled MiXCR 4.7.0
example. Exported column names depend on the MiXCR version and export
options.

## Usage

<pre><code class='language-R'>read_mixcr(
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
One or more repertoire file paths or a file pattern such as
<code>“data/\*.tsv”</code>. Use <code>“\<manifest\>”</code> to read the
paths from <code>manifest</code>.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="schema">schema</code>
</td>
<td>
Column names that define receptor identity, such as <code>c(“v_call”,
“junction_aa”)</code>. Use
<code>immundata::make_receptor_schema()</code> to also select one chain
or a pair of chains from each cell.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="manifest">manifest</code>
</td>
<td>
An optional data frame with file paths and sample annotations.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="manifest_file_col">manifest_file_col</code>
</td>
<td>
The name of the manifest column with file paths when <code>path =
“\<manifest\>”</code>.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="column_profile">column_profile</code>
</td>
<td>
Which input annotations to keep. <code>“minimal”</code> keeps core
identifiers, gene calls, junctions, and counts. <code>“standard”</code>
also keeps sequence, region, and alignment annotations.
<code>“full”</code> keeps all input annotations after the changes
described under Value. Columns required by <code>schema</code> and
<code>repertoire_schema</code> are kept with every profile.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="repertoire_schema">repertoire_schema</code>
</td>
<td>
Columns that define repertoire groups. This is passed to
<code>immundata::read_repertoires()</code>.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="...">…</code>
</td>
<td>
Additional arguments passed to
<code>immundata::read_repertoires()</code>.
</td>
</tr>
</table>

## Value

An immundata::ImmunData object containing the imported MiXCR single-cell
clones. The reader copies source fields to the AIRR output column names
below when those output columns are not already present. Existing AIRR
output column names are kept. In the full profile, source columns remain
alongside the output columns.

<code>cellGroup</code> is copied to <code>cell_id</code>, and
<code>topChains</code> is copied to <code>locus</code>.
<code>allVHitsWithScore</code>, <code>allDHitsWithScore</code>, and
<code>allJHitsWithScore</code> are copied to <code>v_call</code>,
<code>d_call</code>, and <code>j_call</code>. Scores in parentheses,
such as <code>(123)</code>, are removed from the output gene values,
including values in existing gene output columns. <code>nSeqCDR3</code>
is copied to <code>junction</code>, and <code>aaSeqCDR3</code> to
<code>junction_aa</code>. <code>uniqueMoleculeCount</code> is copied to
<code>umi_count</code>, and <code>readCount</code> to
<code>consensus_count</code>. The reader does not copy
<code>readCount</code> to <code>duplicate_count</code> or infer that
count from reads. <strong>Column names kept unchanged:</strong>
<code>duplicate_count</code> and <code>productive</code> keep their
names when present. Rows that are not productive are filtered out when
that column is supplied. Other source columns keep their names in the
full profile.

## Examples

``` r
library("immunread")

# Use one thread for this small example.
duckplyr::db_exec("SET threads TO 1")

path <- system.file("extdata", "mixcr", "4.7.0", "contigs.csv.gz",
  package = "immunread", mustWork = TRUE
)
idata <- read_mixcr(path,
  schema = c("v_call", "junction_aa"),
  output_folder = tempfile("immunread-")
)
```
