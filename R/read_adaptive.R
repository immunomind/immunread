#' @title Read Adaptive immunoSEQ Analyzer and Connector tables
#' @description Reads Adaptive immunoSEQ Analyzer or Connector tables and maps
#'   their resolved gene, CDR3, sequence, and abundance fields to AIRR-style columns.
#'
#' Resolved gene calls are used directly. `templates` supplies bulk abundance
#' for Analyzer and Connector; Analyzer v2 uses `count (templates/reads)`.
#' `seq_reads` remains a distinct source count when present.
#' @inheritParams read_airr_bulk
#' @examples
#' path <- system.file("extdata", "adaptive", "analyzer", "contigs.tsv.gz",
#'   package = "immunread", mustWork = TRUE)
#' idata <- read_adaptive(path, schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-"))
#' @export
read_adaptive <- function(path, schema, manifest = NULL,
                          manifest_file_col = "file", column_profile = "standard",
                          repertoire_schema = "<auto>", ...) {
  prepare <- function(data, manifest = NULL) {
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
      data <- format_aliases(data, c(
        v_call = "v_resolved", d_call = "d_resolved", j_call = "j_resolved",
        junction = junction_col, junction_aa = "cdr3_amino_acid",
        duplicate_count = "templates"
      ))
    }
    if ("frame_type" %in% names(data)) {
      data <- dplyr::mutate(data, productive = frame_type == "In")
    }
    profile_columns(
      keep_productive(data), schema, column_profile,
      repertoire_schema
    )
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    count_col = "duplicate_count", preprocess = list(adaptive = prepare),
    rename_columns = NULL, repertoire_schema = repertoire_schema, ...
  )
}
