

# Read 10x Genomics V(D)J contigs

[**Source code**](https://github.com/immunomind/immunread/tree/main/R/read_10x.R#L45)

## Description

Reads 10x Genomics V(D)J
<a href="https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-5p-outputs-annotations-vdj">contig
annotation CSV files</a> such as
<code>filtered_contig_annotations.csv</code> and maps their barcode,
locus, gene, sequence, and count fields to AIRR-style columns.

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

An immundata::ImmunData object containing the imported 10x contigs.
Source columns are copied to the AIRR output column names below when
those output columns are not already present. Existing AIRR output
column names are kept. The full profile also includes source columns,
except where the region and CDR3 fields are renamed or removed as
described below.

The reader maps <code>barcode</code> to <code>cell_id</code>,
<code>contig_id</code> to <code>sequence_id</code>, and
<code>chain</code> to <code>locus</code>. It maps <code>v_gene</code>,
<code>d_gene</code>, <code>j_gene</code>, and <code>c_gene</code> to
<code>v_call</code>, <code>d_call</code>, <code>j_call</code>, and
<code>c_call</code>. It maps <code>cdr3_nt</code> to
<code>junction</code>, <code>cdr3</code> to <code>junction_aa</code>,
<code>umis</code> to <code>umi_count</code>, and <code>reads</code> to
<code>consensus_count</code>.

In 10x files, <code>fwr1</code> contains amino acid sequences, while
<code>fwr1_nt</code> contains nucleotide sequences. The reader copies
<code>fwr1</code> to <code>fwr1_aa</code>, then replaces
<code>fwr1</code> with the values from <code>fwr1_nt</code>. If
<code>fwr1_nt</code> is absent, <code>fwr1</code> is dropped after its
values are copied to <code>fwr1_aa</code>. The same rule applies to
<code>cdr1</code>, <code>fwr2</code>, <code>cdr2</code>,
<code>fwr3</code>, and <code>fwr4</code>. These copies overwrite any
existing values in the destination columns. The CDR3 source columns are
renamed <code>source_cdr3</code> and <code>source_cdr3_nt</code> when
present. If <code>productive</code> is supplied, unproductive contigs
are removed. <code>umi_count</code> records distinct UMIs, and
<code>consensus_count</code> records supporting reads. Other source
columns keep their names in the full profile. <code>immundata</code>
selects chains according to <code>schema</code>, using
<code>umi_count</code> to rank them.

## Examples

``` r
library("immunread")

# Use one thread for this small example.
duckplyr::db_exec("SET threads TO 1")

path <- system.file("extdata", "10x", "9.0.1", "contigs.csv.gz",
  package = "immunread", mustWork = TRUE
)
idata <- read_10x(path,
  schema = c("v_call", "junction_aa"),
  output_folder = tempfile("immunread-")
)
```
