#' @title Read MiXCR single-cell clones
#' @description Reads MiXCR
#'   [single-cell clone exports](https://mixcr.com/mixcr/reference/mixcr-export/#single-cell-libraries)
#'   and maps their cell, chain, gene,
#'   sequence, and count fields to AIRR-style columns.
#'
#' The input must have `cellGroup`, `topChains`, and `uniqueMoleculeCount`
#' columns, as in the bundled MiXCR 4.7.0 example. Exported column names depend on
#' the MiXCR version and export options.
#'
#' @inheritParams read_airr_sc
#' @return An [immundata::ImmunData] object containing the imported MiXCR
#'   single-cell clones. The reader copies source fields to the AIRR output
#'   column names below when those output columns are not already present.
#'   Existing AIRR output column names are kept. In the full profile, source
#'   columns remain alongside the output columns.
#'
#'   `cellGroup` is copied to `cell_id`, and `topChains` is copied to `locus`.
#'   `allVHitsWithScore`, `allDHitsWithScore`, and `allJHitsWithScore` are copied
#'   to `v_call`, `d_call`, and `j_call`. Scores in parentheses, such as `(123)`,
#'   are removed from the output gene values, including values in existing gene output
#'   columns. `nSeqCDR3` is copied to `junction`, and `aaSeqCDR3` to
#'   `junction_aa`. `uniqueMoleculeCount` is copied to `umi_count`, and
#'   `readCount` to `consensus_count`. The reader does not copy `readCount` to
#'   `duplicate_count` or infer that count from reads. **Column names kept
#'   unchanged:** `duplicate_count` and `productive` keep their names when
#'   present. Rows that are not productive are filtered out when that column is
#'   supplied. Other source columns keep their names in the full profile.
#' @examples
#' # Use one thread for this small example.
#' duckplyr::db_exec("SET threads TO 1")
#'
#' path <- system.file("extdata", "mixcr", "4.7.0", "contigs.csv.gz",
#'   package = "immunread", mustWork = TRUE
#' )
#' idata <- read_mixcr(path,
#'   schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-")
#' )
#' @export
read_mixcr <- function(path, schema, manifest = NULL,
                       manifest_file_col = "file", column_profile = "standard",
                       repertoire_schema = "<auto>", ...) {
  prepare <- function(data, manifest = NULL) {
    data <- format_aliases(data, c(
      cell_id = "cellGroup", locus = "topChains",
      v_call = "allVHitsWithScore", d_call = "allDHitsWithScore",
      j_call = "allJHitsWithScore", junction = "nSeqCDR3",
      junction_aa = "aaSeqCDR3", umi_count = "uniqueMoleculeCount",
      consensus_count = "readCount"
    ))
    for (name in intersect(c("v_call", "d_call", "j_call"), names(data))) {
      data <- dplyr::mutate(
        data,
        !!name := gsub("\\([^)]*\\)", "", !!rlang::sym(name))
      )
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
    preprocess = list(mixcr = prepare), rename_columns = NULL,
    repertoire_schema = repertoire_schema, ...
  )
}
