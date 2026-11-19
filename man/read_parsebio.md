

# Read Parse Biosciences TCR tables

[**Source code**](https://github.com/immunomind/immunread/tree/main/R/read_parse.R#L53)

## Description

Reads Parse Biosciences TCR tables in which each row describes one cell.
The reader separates the chains into rows and uses AIRR names for the
output annotations.

See Parse's
<a href="https://resources.parsebiosciences.com/support-suite/hc/en-us/articles/44260645430676-understanding-the-anndata-scanpy-object-that-is-downloadable-from-trailmaker-insights-module">public
immune annotation overview</a> for descriptions of the wide per-cell
fields. This reader expects the wide layout with <code>Barcode</code>,
<code style="white-space: pre;">TRA\_*</code>,
<code style="white-space: pre;">TRB\_*</code>, and optional
<code style="white-space: pre;">secondary_TRA\_*</code> /
<code style="white-space: pre;">secondary_TRB\_*</code> columns. For
AIRR-formatted <code>tcr_annotation_airr.tsv</code> exports, use
<code>read_airr_sc()</code>.

The reader prepares one row for each supplied primary or secondary TRA
or TRB chain. Rows without a CDR3 amino acid sequence are excluded
before 'immundata' selects the chains used to construct receptors.

## Usage

<pre><code class='language-R'>read_parsebio(
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
available chain annotations. For each of <code>TRA</code>,
<code>TRB</code>, <code>secondary_TRA</code>, and
<code>secondary_TRB</code>, the reader copies
<code style="white-space: pre;">\*\_V</code> to <code>v_call</code>,
<code style="white-space: pre;">\*\_D</code> to <code>d_call</code>,
<code style="white-space: pre;">\*\_J</code> to <code>j_call</code>, and
<code style="white-space: pre;">\*\_C</code> to <code>c_call</code>. It
copies <code style="white-space: pre;">\*\_cdr3_aa</code> to
<code>junction_aa</code> and
<code style="white-space: pre;">\*\_transcript_count</code> to
<code>umi_count</code> for chain ranking and
<code style="white-space: pre;">\*\_read_count</code> to
<code>consensus_count</code>. <code>Barcode</code> maps to
<code>cell_id</code>. <code>locus</code> is set to <code>TRA</code> for
primary and secondary TRA chains, and to <code>TRB</code> for primary
and secondary TRB chains. The group name is recorded in
<code>source_chain</code>.

AIRR output columns that are already present keep their names and are
used instead of the corresponding source columns. The reader always sets
<code>locus</code> and <code>source_chain</code> from the chain group.

If no group has a read-count column, <code>consensus_count</code> is
absent. If some groups have that column, chains from groups without it
have <code>NA</code> for <code>consensus_count</code>. Missing optional
gene calls also have <code>NA</code>.

With <code>column_profile = "full"</code>, the reader retains the
original columns, including <code>Barcode</code> and the wide chain
columns, and adds the AIRR columns. With <code>column_profile =
"minimal"</code> or <code>"standard"</code>, the reader keeps the
columns included in that profile. The available columns also depend on
the input and <code>schema</code>.

## Examples

``` r
library("immunread")

# Use one thread for this small example.
duckplyr::db_exec("SET threads TO 1")

path <- system.file("extdata", "parsebio", "1.6.1", "contigs.csv.gz",
  package = "immunread", mustWork = TRUE
)
idata <- read_parsebio(path,
  schema = c("v_call", "junction_aa"),
  output_folder = tempfile("immunread-")
)
```
