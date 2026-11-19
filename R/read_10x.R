#' @title Read 10x Genomics V(D)J contigs
#' @description Reads 10x Genomics V(D)J
#'   [contig annotation CSV files](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-5p-outputs-annotations-vdj)
#'   such as `filtered_contig_annotations.csv` and maps their
#'   barcode, locus, gene, sequence, and count fields to AIRR-style columns.
#'
#' @inheritParams read_airr_sc
#' @return An [immundata::ImmunData] object containing the imported 10x
#'   contigs. Source columns are copied to the AIRR output column names below
#'   when those output columns are not already present. Existing AIRR output
#'   column names are kept. The full profile also includes source columns,
#'   except where the region and CDR3 fields are renamed or removed as described
#'   below.
#'
#'   The reader maps `barcode` to `cell_id`, `contig_id` to `sequence_id`, and
#'   `chain` to `locus`. It maps `v_gene`, `d_gene`, `j_gene`, and `c_gene` to
#'   `v_call`, `d_call`, `j_call`, and `c_call`. It maps `cdr3_nt` to
#'   `junction`, `cdr3` to `junction_aa`, `umis` to `umi_count`, and `reads` to
#'   `consensus_count`.
#'
#'   In 10x files, `fwr1` contains amino acid sequences, while `fwr1_nt`
#'   contains nucleotide sequences. The reader copies `fwr1` to `fwr1_aa`,
#'   then replaces `fwr1` with the values from `fwr1_nt`. If `fwr1_nt` is
#'   absent, `fwr1` is dropped after its values are copied to `fwr1_aa`.
#'   The same rule applies to `cdr1`, `fwr2`, `cdr2`, `fwr3`, and `fwr4`.
#'   These copies overwrite any existing values in the destination columns.
#'   The CDR3 source columns are renamed `source_cdr3` and `source_cdr3_nt`
#'   when present. If `productive` is supplied, unproductive contigs are removed.
#'   `umi_count` records distinct UMIs, and `consensus_count` records supporting
#'   reads. Other source columns keep their names in the full profile.
#'   `immundata` selects chains according to `schema`, using `umi_count`
#'   to rank them.
#' @examples
#' # Use one thread for this small example.
#' duckplyr::db_exec("SET threads TO 1")
#'
#' path <- system.file("extdata", "10x", "9.0.1", "contigs.csv.gz",
#'   package = "immunread", mustWork = TRUE
#' )
#' idata <- read_10x(path,
#'   schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-")
#' )
#' @export
read_10x <- function(path, schema, manifest = NULL,
                     manifest_file_col = "file", column_profile = "standard",
                     repertoire_schema = "<auto>", ...) {
  prepare <- function(data, manifest = NULL) {
    data <- format_aliases(data, c(
      cell_id = "barcode", sequence_id = "contig_id", locus = "chain",
      v_call = "v_gene", d_call = "d_gene", j_call = "j_gene",
      c_call = "c_gene", junction = "cdr3_nt",
      junction_aa = "cdr3", umi_count = "umis", consensus_count = "reads"
    ))
    # 10x's unsuffixed region fields are amino acids; AIRR's are nucleotides.
    # Very incovenient, but what can we do...
    for (region in c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "fwr4")) {
      if (!region %in% names(data)) next
      data <- dplyr::mutate(data, !!paste0(region, "_aa") := !!rlang::sym(region))
      nt <- paste0(region, "_nt")
      if (nt %in% names(data)) {
        data <- dplyr::mutate(data, !!region := !!rlang::sym(nt))
      } else {
        data <- dplyr::select(data, -dplyr::all_of(region))
      }
    }
    if ("cdr3" %in% names(data)) data <- dplyr::rename(data, source_cdr3 = cdr3)
    if ("cdr3_nt" %in% names(data)) {
      data <- dplyr::rename(data, source_cdr3_nt = cdr3_nt)
    }
    profile_columns(
      keep_productive(data), schema, column_profile,
      repertoire_schema,
      count_col = "umi_count"
    )
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = "cell_id", locus_col = "locus", umi_col = "umi_count",
    preprocess = list(`10x` = prepare), rename_columns = NULL,
    repertoire_schema = repertoire_schema, ...
  )
}
