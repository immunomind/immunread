

# Read Adaptive immunoSEQ Analyzer and Connector tables

[**Source code**](https://github.com/immunomind/immunread/tree/main/R/read_adaptive.R#L51)

## Description

Reads Adaptive immunoSEQ Analyzer or Connector tables.

Input columns are described in Adaptive's
<a href="https://www.adaptivebiotech.com/wp-content/uploads/2019/07/MRK-00342_immunoSEQ_TechNote_DataExport_WEB_REV.pdf">Analyzer
TSV export guide</a> and the legacy
<a href="https://clients.adaptivebiotech.com/assets/downloads/immunoSEQ_AnalyzerManual.pdf">Analyzer
manual, Sample Export section</a>.

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

An immundata::ImmunData object containing the imported Adaptive Analyzer
or Connector rearrangements. The reader copies source fields to the AIRR
output column names below when those output columns are not already
present. Existing AIRR output column names are kept. The full profile
also includes source columns, except that
<code>rearrangement_type</code> is renamed to
<code>source_rearrangement_type</code>.

In Analyzer v2 files identified by <code>sequenceStatus</code>, the
reader copies <code>vMaxResolved</code>, <code>dMaxResolved</code>, and
<code>jMaxResolved</code> to <code>v_call</code>, <code>d_call</code>,
and <code>j_call</code>. It copies <code>aminoAcid</code> to
<code>junction_aa</code>, <code>nucleotide</code> to
<code>sequence</code>, and <code>count (templates/reads)</code> to
<code>duplicate_count</code>. The reader sets <code>productive</code>
from <code>sequenceStatus == "In"</code>. When <code>frame_type</code>
is present, it sets <code>productive</code> from <code>frame_type ==
"In"</code> instead.

In Analyzer and Connector files without <code>sequenceStatus</code>, the
reader copies <code>v_resolved</code>, <code>d_resolved</code>, and
<code>j_resolved</code> to <code>v_call</code>, <code>d_call</code>, and
<code>j_call</code>. It copies the first available column among
<code>cdr3_rearrangement</code> and <code>cdr3_nucleotide</code> to
<code>junction</code>, and copies <code>cdr3_amino_acid</code> to
<code>junction_aa</code>. It copies the first available column among
<code>rearrangement</code> and <code>nucleotide</code> to
<code>sequence</code>, and copies <code>templates</code> to
<code>duplicate_count</code>.

Nonproductive rows are removed. Source <code>seq_reads</code>,
<code>sequenceStatus</code>, and <code>frame_type</code> keep their
names in the full profile, as do other source fields. The source
<code>rearrangement_type</code> describes gene composition and is
retained as <code>source_rearrangement_type</code> in the full profile.
This is not AIRR's field for observed, simulated, or inferred origin.
Adaptive counts represent copies, not distinct UMIs. The reader does not
calculate UMI counts.

## Examples

``` r
library("immunread")

# Use one thread for this small example.
duckplyr::db_exec("SET threads TO 1")

path <- system.file("extdata", "adaptive", "analyzer", "contigs.tsv.gz",
  package = "immunread", mustWork = TRUE
)
idata <- read_adaptive(path,
  schema = c("v_call", "junction_aa"),
  output_folder = tempfile("immunread-")
)
```
