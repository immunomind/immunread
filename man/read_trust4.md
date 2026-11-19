

# Read TRUST4 cell tables

[**Source code**](https://github.com/immunomind/immunread/tree/main/R/read_trust4.R#L49)

## Description

Reads TRUST4 cell tables and separates their chain records into rows
with AIRR names for the gene, sequence, and count annotations.

The
<a href="https://github.com/liulab-dfci/TRUST4#10x-genomics-data-and-barcode-based-single-cell-data">TRUST4
barcode report format</a> describes the embedded comma-separated chain
fields. This reader expects the tabular layout in the bundled example,
with <code>X.barcode</code>, <code>chain1</code>, <code>chain2</code>,
<code>secondary_chain1</code>, and <code>secondary_chain2</code>
columns. See the TRUST4 paper,
[doi:10.1038/s41592-021-01142-2](https://doi.org/10.1038/s41592-021-01142-2).

Each non-empty primary or secondary chain record becomes a row before
'immundata' selects the chains used to construct receptors.

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
id="...">...</code></td>
<td>Additional arguments passed to
<code>immundata::read_repertoires()</code>.</td>
</tr>
</tbody>
</table>

## Value

An immundata::ImmunData object containing the imported repertoires and
available chain annotations. The reader splits <code>chain1</code>,
<code>chain2</code>, <code>secondary_chain1</code>, and
<code>secondary_chain2</code> at semicolons and makes each record into
its own row. Within each record, comma-separated fields 1 through 8 map
in order to <code>v_call</code>, <code>d_call</code>,
<code>j_call</code>, <code>c_call</code>, <code>junction</code>,
<code>junction_aa</code>, <code>source_read_count</code>, and
<code>sequence_id</code>. <code>X.barcode</code> maps to
<code>cell_id</code>, and the field containing the record is recorded as
<code>source_chain</code>. The <code>locus</code> value comes from the
first three characters of <code>v_call</code>, or from
<code>j_call</code> when <code>v_call</code> is <code>\*</code>.

The read-support estimate in field 7 becomes
<code>source_read_count</code>. Its fractional part is removed to give
<code>consensus_count</code>, matching the upstream AIRR exporter. For
example, a read-support estimate of 7.9 gives
<code>consensus_count</code> 7. TRUST4 provides no separate UMI count,
so <code>consensus_count</code> is used for chain ranking. Empty and
<code>\*</code> chain fields are skipped. With <code>column_profile =
"full"</code>, the reader retains the original columns, including
<code>X.barcode</code> and the embedded chain fields, and adds the AIRR
columns. <code>source_read_count</code> is also available with this
profile. The contents of the chain fields are split into the mapped
columns above. With <code>column_profile = "minimal"</code> or
<code>"standard"</code>, the reader keeps the columns included in that
profile. The available columns also depend on the input and
<code>schema</code>.

## Examples

``` r
library("immunread")

# Use one thread for this small example.
duckplyr::db_exec("SET threads TO 1")

path <- system.file("extdata", "trust4", "1.1.7", "contigs.csv.gz",
  package = "immunread", mustWork = TRUE
)
idata <- read_trust4(path,
  schema = c("v_call", "junction_aa"),
  output_folder = tempfile("immunread-")
)
```
