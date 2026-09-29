#' @title Read BD Rhapsody AIRR-like contigs
#' @description Reads BD Rhapsody AIRR-like contigs and preserves their cell,
#'   locus, receptor, and count fields in an `ImmunData` object.
#'
#' Supplied secondary chains are retained; `immundata` ranks chains when
#' constructing receptors. BD count fields keep their source meanings.
#' @inheritParams read_airr_sc
#' @examples
#' path <- system.file("extdata", "bd", "2.3", "contigs.csv.gz",
#'   package = "immunread", mustWork = TRUE)
#' idata <- read_bd(path, schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-"))
#' @export
read_bd <- function(path, schema, manifest = NULL,
                    manifest_file_col = "file", column_profile = "standard",
                    repertoire_schema = "<auto>", ...) {
  prepare <- function(data, manifest = NULL) {
    profile_columns(
      keep_productive(data), schema, column_profile,
      repertoire_schema
    )
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = "cell_id", locus_col = "locus", umi_col = "duplicate_count",
    preprocess = list(profile = prepare), rename_columns = NULL,
    repertoire_schema = repertoire_schema, ...
  )
}
