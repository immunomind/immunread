

# Read an AIRR Rearrangement table as single-cell repertoire data

[**Source code**](https://github.com/immunomind/immunread/tree/main/R/read_airr.R#L143)

## Description

Reads a single-cell
<a href="https://docs.airr-community.org/en/stable/datarep/rearrangements.html">AIRR
Rearrangement table</a> and uses its cell identifiers and counts to
select receptor chains.

To select chains, the function uses <code>umi_count</code> when it is
present. If that column is absent, it uses <code>duplicate_count</code>.
If both are absent, it uses <code>consensus_count</code>. This is the
same UMI-first priority used by <code>read_airr_bulk()</code>, but here
the counts are used to rank chains within each cell. Missing values in
the selected field cause an error even when another count field has
values. The count and missing-value rules described in
<code>read_airr_bulk()</code> also apply here.

## Usage

<pre><code class='language-R'>read_airr_sc(
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

An immundata::ImmunData object containing the imported single-cell
repertoire data. If <code>cell_id</code> is absent, values from
<code>barcode</code> are copied into a new <code>cell_id</code> column.
If <code>cell_id</code> is already present, it is used without changing
its name. When a source <code>barcode</code> column is supplied, it is
also available with <code>column_profile = "full"</code>.

Other AIRR columns keep their original names. This includes gene calls,
<code>locus</code>, <code>junction</code>, <code>junction_aa</code>,
<code>duplicate_count</code>, <code>umi_count</code>, and
<code>consensus_count</code> when present. The function does not rename
a count field when it uses that field to select chains. These count
names describe sequence copies, distinct UMIs, and reads used to build a
consensus or contig in AIRR terminology. The reader also keeps
<code>cdr3</code> and <code>cdr3_aa</code> separate from the junction
fields.

The available annotations depend on the input and
<code>column_profile</code>. 'immundata' adds its own analysis columns,
and <code>manifest</code> can add sample annotations under their
original column names.

## Examples

``` r
library("immunread")

# Use one thread for this small example.
duckplyr::db_exec("SET threads TO 1")

path <- system.file("extdata", "airr", "2.0-sc", "contigs.tsv.gz",
  package = "immunread", mustWork = TRUE
)
schema <- immundata::make_receptor_schema(
  c("v_call", "junction_aa"),
  chains = "TRB"
)
idata <- read_airr_sc(path,
  schema = schema,
  output_folder = tempfile("immunread-")
)
```
