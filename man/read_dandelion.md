

# Read Dandelion AIRR-like contigs

[**Source code**](https://github.com/immunomind/immunread/tree/main/R/read_dandelion.R#L36)

## Description

Reads Dandelion
<a href="https://sc-dandelion.readthedocs.io/en/latest/notebooks/base/0_dandelion_primer.html">AIRR
contig tables</a> (the <code>Dandelion.data</code> structure) and reads
their cell, locus, receptor, and abundance fields into an
immundata::ImmunData object.

The reader passes available chains to <code>immundata</code>, which
selects chains according to <code>schema</code> and count ranking. It
ranks chains by <code>umi_count</code> when present. If that field is
absent, it uses <code>duplicate_count</code>. If both fields are absent,
it uses <code>consensus_count</code>. Counts used for ranking must be
finite whole numbers that are zero or greater, with no missing values.

## Usage

<pre><code class='language-R'>read_dandelion(
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

An immundata::ImmunData object containing the imported repertoires and
available annotations. <strong>Column names kept unchanged:</strong>
Supplied AIRR names such as <code>cell_id</code>, <code>locus</code>,
<code>v_call</code>, <code>d_call</code>, <code>j_call</code>,
<code>junction</code>, <code>junction_aa</code>, and
<code>consensus_count</code> keep their names. When
<code>umi_count</code> is absent and <code>duplicate_count</code> is
present, the source <code>duplicate_count</code> column is renamed to
<code>umi_count</code>. Its original value is also kept as
<code>source_duplicate_count</code> in the <code>full</code> column
profile. When <code>umi_count</code> is already present,
<code>duplicate_count</code> keeps its name and its separate copy count.
<code>umi_count</code> counts distinct UMIs and
<code>consensus_count</code> counts reads supporting consensus
sequences. Other available columns depend on the input and
<code>column_profile</code>. Source and manifest annotations keep their
names.

## Examples

``` r
library("immunread")

path <- system.file("extdata", "dandelion", "0.5.5", "contigs.csv.gz",
  package = "immunread", mustWork = TRUE
)
idata <- read_dandelion(path,
  schema = c("v_call", "junction_aa"),
  output_folder = tempfile("immunread-")
)
```
