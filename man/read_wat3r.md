

# Read WAT3R TCR tables

[**Source code**](https://github.com/immunomind/immunread/tree/main/R/read_wat3r.R#L56)

## Description

Reads WAT3R
<a href="https://github.com/mainciburu/WAT3R#downstream">barcode_results.csv
tables</a> in which each row describes one cell. The reader separates
the chains into rows and uses AIRR names for the output annotations.

The input contains <code>BC</code>,
<code style="white-space: pre;">TRA\_*</code>,
<code style="white-space: pre;">TRB\_*</code>, and optional
<code style="white-space: pre;">TRA.2\_\*</code> columns. See the WAT3R
paper,
[doi:10.1093/bioinformatics/btac382](https://doi.org/10.1093/bioinformatics/btac382).

The reader prepares one row for each supplied TRA, TRB, or secondary TRA
chain. Rows without a CDR3 amino acid sequence are excluded before
‘immundata’ selects the chains used to construct receptors.

## Usage

<pre><code class='language-R'>read_wat3r(
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
available chain annotations. Gene column names are converted for each
chain group as follows:

<ul>
<li>

For TRA, <code>TRAV</code> becomes <code>v_call</code> and
<code>TRAJ</code> becomes <code>j_call</code>.

</li>
<li>

For TRB, <code>TRBV</code> becomes <code>v_call</code>,
<code>TRBD</code> becomes <code>d_call</code>, and <code>TRBJ</code>
becomes <code>j_call</code>.

</li>
<li>

For secondary TRA, <code>TRAV.2</code> becomes <code>v_call</code> and
<code>TRAJ.2</code> becomes <code>j_call</code>.

</li>
</ul>

The matching <code>TRA_CDR3nuc</code>, <code>TRB_CDR3nuc</code>, and
<code>TRA.2_CDR3nuc</code> values map to <code>cdr3</code>. These
nucleotide sequences exclude conserved anchor codons. The matching
<code style="white-space: pre;">\*\_CDR3</code> amino acid sequences map
to <code>junction_aa</code> and include the conserved amino acids. WAT3R
does not supply a nucleotide <code>junction</code>, so this reader does
not create one. The matching
<code style="white-space: pre;">\*\_nReads</code> values map to
<code>consensus_count</code>.
<code style="white-space: pre;">\*\_CDR3_UMIcount</code> values map to
<code>umi_count</code> for chain ranking. <code>BC</code> maps to
<code>cell_id</code>. The chain group supplies <code>locus</code>
(<code>TRA.2</code> is recorded as <code>TRA</code>) and
<code>source_chain</code>.

AIRR output columns that are already present keep their names and are
used instead of the corresponding source columns. The reader always sets
<code>locus</code> and <code>source_chain</code> from the chain group.

If no group has an <code>nReads</code> column,
<code>consensus_count</code> is absent. If some groups have that column,
chains from groups without it have <code>NA</code> for
<code>consensus_count</code>. Missing optional gene calls also have
<code>NA</code>.

With <code>column_profile = “full”</code>, original source columns,
including <code>BC</code> and the wide chain columns, remain under their
original names and the mapped columns are added. <strong>These source
column names are kept unchanged.</strong> With <code>column_profile =
“minimal”</code> or <code>“standard”</code>, the reader keeps the
columns included in that profile. The available columns also depend on
the input and <code>schema</code>.

## Examples

``` r
library("immunread")

# Use one thread for this small example.
duckplyr::db_exec("SET threads TO 1")

path <- system.file("extdata", "wat3r", "1.3", "contigs.csv.gz",
  package = "immunread", mustWork = TRUE)
idata <- read_wat3r(path, schema = c("v_call", "junction_aa"),
  output_folder = tempfile("immunread-"))
```
