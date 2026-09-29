#' @title Read Dandelion AIRR-like contigs
#' @description Reads Dandelion AIRR-like contigs and maps their cell, locus,
#'   receptor, and abundance fields to AIRR-style columns.
#'
#' Supplied secondary chains are retained; `immundata` ranks chains when
#' constructing receptors.
#' @inheritParams read_airr_sc
#' @examples
#' path <- system.file("extdata", "dandelion", "0.5.5", "contigs.csv.gz",
#'   package = "immunread", mustWork = TRUE)
#' idata <- read_dandelion(path, schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-"))
#' @export
read_dandelion <- function(path, schema, manifest = NULL,
                           manifest_file_col = "file", column_profile = "standard",
                           repertoire_schema = "<auto>", ...) {
  columns <- source_columns(path, manifest, manifest_file_col)
  rank_col <- first_column(
    columns,
    c("umi_count", "duplicate_count", "consensus_count")
  )
  prepare <- function(data, manifest = NULL) {
    profile_columns(
      keep_productive(data), schema, column_profile,
      repertoire_schema
    )
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = "cell_id",
    locus_col = if ("locus" %in% columns) "locus" else NULL,
    umi_col = rank_col, preprocess = list(profile = prepare),
    rename_columns = NULL, repertoire_schema = repertoire_schema, ...
  )
}
