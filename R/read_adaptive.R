#' @title Read Adaptive immunoSEQ Analyzer and Connector tables
#' @description Reads Adaptive immunoSEQ Analyzer or Connector tables.
#'
#' Input columns are described in Adaptive's
#' [Analyzer TSV export guide](https://www.adaptivebiotech.com/wp-content/uploads/2019/07/MRK-00342_immunoSEQ_TechNote_DataExport_WEB_REV.pdf)
#' and the legacy
#' [Analyzer manual, Sample Export section](https://clients.adaptivebiotech.com/assets/downloads/immunoSEQ_AnalyzerManual.pdf).
#'
#' @inheritParams read_airr_bulk
#' @return An [immundata::ImmunData] object containing the imported Adaptive
#'   Analyzer or Connector rearrangements. The reader copies source fields to
#'   the AIRR output column names below when those output columns are not already
#'   present. Existing AIRR output column names are kept. The full profile also
#'   includes source columns, except that `rearrangement_type` is renamed to
#'   `source_rearrangement_type`.
#'
#'   In Analyzer v2 files identified by `sequenceStatus`, the reader copies
#'   `vMaxResolved`, `dMaxResolved`, and `jMaxResolved` to `v_call`, `d_call`,
#'   and `j_call`. It copies `aminoAcid` to `junction_aa`, `nucleotide` to
#'   `sequence`, and `count (templates/reads)` to `duplicate_count`. The reader
#'   sets `productive` from `sequenceStatus == "In"`. When `frame_type` is
#'   present, it sets `productive` from `frame_type == "In"` instead.
#'
#'   In Analyzer and Connector files without `sequenceStatus`, the reader copies
#'   `v_resolved`, `d_resolved`, and `j_resolved` to `v_call`, `d_call`, and
#'   `j_call`. It copies the first available column among
#'   `cdr3_rearrangement` and `cdr3_nucleotide` to `junction`, and copies
#'   `cdr3_amino_acid` to `junction_aa`. It copies the first available column
#'   among `rearrangement` and `nucleotide` to `sequence`, and copies
#'   `templates` to `duplicate_count`.
#'
#'   Nonproductive rows are removed. Source `seq_reads`,
#'   `sequenceStatus`, and `frame_type` keep their names in the full profile,
#'   as do other source fields. The source `rearrangement_type` describes gene
#'   composition and is retained as
#'   `source_rearrangement_type` in the full profile. This is not AIRR's field
#'   for observed, simulated, or inferred origin. Adaptive counts represent
#'   copies, not distinct UMIs. The reader does not calculate UMI counts.
#' @examples
#' # Use one thread for this small example.
#' duckplyr::db_exec("SET threads TO 1")
#'
#' path <- system.file("extdata", "adaptive", "analyzer", "contigs.tsv.gz",
#'   package = "immunread", mustWork = TRUE
#' )
#' idata <- read_adaptive(path,
#'   schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-")
#' )
#' @export
read_adaptive <- function(path, schema, manifest = NULL,
                          manifest_file_col = "file", column_profile = "standard",
                          repertoire_schema = "<auto>", ...) {
  prepare <- function(data, manifest = NULL) {
    if ("rearrangement_type" %in% names(data)) {
      data <- dplyr::rename(data, source_rearrangement_type = rearrangement_type)
    }
    if ("sequenceStatus" %in% names(data)) {
      data <- format_aliases(data, c(
        v_call = "vMaxResolved", d_call = "dMaxResolved",
        j_call = "jMaxResolved", junction_aa = "aminoAcid",
        sequence = "nucleotide",
        duplicate_count = "count (templates/reads)"
      ))
      data <- dplyr::mutate(data, productive = sequenceStatus == "In")
    } else {
      junction_col <- first_column(names(data), c("cdr3_rearrangement", "cdr3_nucleotide"))
      sequence_col <- first_column(names(data), c("rearrangement", "nucleotide"))
      data <- format_aliases(data, c(
        v_call = "v_resolved", d_call = "d_resolved", j_call = "j_resolved",
        junction = junction_col, junction_aa = "cdr3_amino_acid", sequence = sequence_col,
        duplicate_count = "templates"
      ))
    }
    if ("frame_type" %in% names(data)) {
      data <- dplyr::mutate(data, productive = frame_type == "In")
    }
    profile_columns(
      keep_productive(data), schema, column_profile,
      repertoire_schema,
      count_col = "duplicate_count"
    )
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    count_col = "duplicate_count", preprocess = list(adaptive = prepare),
    rename_columns = NULL, repertoire_schema = repertoire_schema, ...
  )
}
