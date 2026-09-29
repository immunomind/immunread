#' @title Read Immcantation AIRR-like tables
#' @description Reads Immcantation AIRR-like tables in bulk or single-cell mode
#'   and maps their receptor and abundance fields to AIRR-style columns.
#'
#' Supplied secondary chains are retained in single-cell mode; `immundata`
#' ranks chains when constructing receptors.
#' @inheritParams read_airr_sc
#' @param mode `"sc"` for single-cell tables or `"bulk"` for bulk tables.
#' @examples
#' path <- system.file("extdata", "immcantation", "4.6.0", "contigs.csv.gz",
#'   package = "immunread", mustWork = TRUE)
#' idata <- read_immcantation(path, schema = c("v_call", "junction_aa"),
#'   mode = "sc", output_folder = tempfile("immunread-"))
#' @export
read_immcantation <- function(path, schema, manifest = NULL,
                              manifest_file_col = "file", column_profile = "standard",
                              repertoire_schema = "<auto>", mode = c("sc", "bulk"),
                              ...) {
  mode <- match.arg(mode)
  columns <- source_columns(path, manifest, manifest_file_col)
  count_col <- if (mode == "bulk") {
    first_column(columns, c("duplicate_count", "umi_count"))
  } else NULL
  rank_col <- if (mode == "sc") {
    first_column(columns, c("umi_count", "duplicate_count", "consensus_count"))
  } else NULL
  prepare <- function(data, manifest = NULL) {
    profile_columns(keep_productive(data), schema, column_profile,
                    repertoire_schema)
  }
  immundata::read_repertoires(
    path, schema, manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = if (mode == "sc") "cell_id" else NULL,
    count_col = count_col,
    locus_col = if ("locus" %in% columns) "locus" else NULL,
    umi_col = rank_col, preprocess = list(profile = prepare),
    rename_columns = NULL, repertoire_schema = repertoire_schema, ...
  )
}
