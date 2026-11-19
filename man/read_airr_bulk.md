

# Read an AIRR Rearrangement table as bulk repertoire data

[**Source code**](https://github.com/immunomind/immunread/tree/main/R/read_airr.R#L78)

## Description

Reads a bulk
<a href="https://docs.airr-community.org/en/stable/datarep/rearrangements.html">AIRR
Rearrangement table</a>.

## Usage

<pre><code class='language-R'>read_airr_bulk(
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

## Details

Counts are taken from the first available column in this order:
<code>umi_count</code>, <code>duplicate_count</code>,
<code>consensus_count</code>. If none is present, no count column is
passed to <code>immundata</code>.

So in short, this default prioritizes counts of tagged molecules over
read counts. Please mind that the
<a href="https://docs.airr-community.org/en/stable/datarep/news.html#version-1-4-1-august-27-2022">AIRR
schema v1.4.1 release notes</a> distinguish distinct UMIs
(<code>umi_count</code>) from duplicate sequence observations
(<code>duplicate_count</code>), and this priority order is our
interpretation of the standard.

Count values must be finite, non-negative integers. Optional count
columns can contain missing values, but the selected count column must
have a value for every productive input row. Missing counts cause an
error; they are not replaced with zero or one.

In recognized sequence and count fields, text such as <code>"NA"</code>
and <code>"None"</code> becomes a missing value (<code>NA</code> in R).
A row is excluded if it has no value for a receptor feature required by
<code>schema</code>. These changes do not apply to other source
annotations or to metadata added from <code>manifest</code>.

## Value

An immundata::ImmunData object containing the imported repertoire data.
The chain annotations keep the supplied AIRR column names, including
<code>v_call</code>, <code>d_call</code>, <code>j_call</code>,
<code>junction</code>, and <code>junction_aa</code> when present.

<code>duplicate_count</code>, <code>umi_count</code>, and
<code>consensus_count</code> remain separate columns with their original
names. In AIRR terminology, they describe sequence copies, distinct
UMIs, and reads used to build a consensus or contig. Selecting a count
column does not rename the other count columns or change their units.
The reader also keeps <code>cdr3</code> and <code>cdr3_aa</code>
separate from the junction fields. It does not use one name as a
replacement for another.

The available annotations depend on the input and
<code>column_profile</code>. 'immundata' adds its own analysis columns,
and <code>manifest</code> can add sample annotations under their
original column names.

## Examples

``` r
library("immunread")

# Use one thread for this small example.
duckplyr::db_exec("SET threads TO 1")

path <- system.file("extdata", "airr", "2.0-bulk", "contigs.tsv.gz",
  package = "immunread", mustWork = TRUE
)
idata <- read_airr_bulk(path,
  schema = c("v_call", "junction_aa"),
  output_folder = tempfile("immunread-")
)
```
