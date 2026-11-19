

# Read Immcantation AIRR-like tables

[**Source code**](https://github.com/immunomind/immunread/tree/main/R/read_immcantation.R#L38)

## Description

Reads Immcantation tables using the
<a href="https://changeo.readthedocs.io/en/stable/standard.html#airr-c-format">AIRR-C
format</a> in bulk or single-cell mode and reads their annotations into
an immundata::ImmunData object.

In single-cell mode, the reader passes available chains to
<code>immundata</code>, which selects chains according to
<code>schema</code> and count ranking. In bulk mode, counts are taken
from <code>duplicate_count</code>, or from <code>umi_count</code> if
<code>duplicate_count</code> is absent. In single-cell mode, chains are
ranked by <code>umi_count</code> when present, then by
<code>duplicate_count</code> if <code>umi_count</code> is absent, and
finally by <code>consensus_count</code> if both are absent. The selected
count column must contain finite, non-negative integers with no missing
values.

## Usage

<pre><code class='language-R'>read_immcantation(
  path,
  schema,
  manifest = NULL,
  manifest_file_col = "file",
  column_profile = "standard",
  repertoire_schema = "&lt;auto&gt;",
  mode = c("sc", "bulk"),
  ...
)
</code></pre>

## Arguments

<table data-quarto-postprocess="true" role="presentation">
<tbody>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="path">path</code></td>
<td>One or more repertoire file paths or a file pattern such as
<code>"data/*.tsv"</code>. Use <code>"&lt;manifest&gt;"</code> to read
the paths from <code>manifest</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="schema">schema</code></td>
<td>Column names that define receptor identity, such as
<code>c("v_call", "junction_aa")</code>. Use
<code>immundata::make_receptor_schema()</code> to also select one chain
or a pair of chains from each cell.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="manifest">manifest</code></td>
<td>An optional data frame with file paths and sample annotations.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="manifest_file_col">manifest_file_col</code></td>
<td>The name of the manifest column with file paths when
<code>path = "&lt;manifest&gt;"</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="column_profile">column_profile</code></td>
<td>Which input annotations to keep. <code>"minimal"</code> keeps core
identifiers, gene calls, junctions, and counts. <code>"standard"</code>
also keeps sequence, region, and alignment annotations.
<code>"full"</code> keeps all input annotations after the changes
described under Value. Columns required by <code>schema</code> and
<code>repertoire_schema</code> are kept with every profile.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="repertoire_schema">repertoire_schema</code></td>
<td>Columns that define repertoire groups. This is passed to
<code>immundata::read_repertoires()</code>.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="mode">mode</code></td>
<td><code>"sc"</code> for single-cell tables or <code>"bulk"</code> for
bulk tables.</td>
</tr>
<tr>
<td
style="white-space: nowrap; font-family: monospace; vertical-align: top"><code
id="...">...</code></td>
<td>Additional arguments passed to
<code>immundata::read_repertoires()</code>.</td>
</tr>
</tbody>
</table>

## Value

An immundata::ImmunData object containing the imported repertoires and
available annotations. Columns retained by <code>column_profile</code>
keep their original names, including <code>cell_id</code>,
<code>locus</code>, <code>v_call</code>, <code>d_call</code>,
<code>j_call</code>, <code>junction</code>, <code>junction_aa</code>,
<code>duplicate_count</code>, <code>umi_count</code>, and
<code>consensus_count</code>. These count columns remain separate.
<code>umi_count</code> counts distinct UMIs,
<code>duplicate_count</code> counts copies, and
<code>consensus_count</code> counts reads supporting consensus
sequences. The reader does not treat <code>duplicate_count</code> as a
UMI count. Other available columns depend on the input and
<code>column_profile</code>. Missing source fields remain absent. Source
and manifest annotations keep their names.

## Examples

``` r
library("immunread")

path <- system.file("extdata", "immcantation", "single-cell.tsv.gz",
  package = "immunread", mustWork = TRUE
)
idata <- read_immcantation(path,
  schema = c("v_call", "junction_aa"),
  mode = "sc", output_folder = tempfile("immunread-")
)
```
