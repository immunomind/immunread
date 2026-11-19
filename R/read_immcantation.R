#' @title Read Immcantation AIRR-C bulk tables
#' @description Reads bulk Immcantation tables using the
#'   [AIRR-C format](https://changeo.readthedocs.io/en/stable/standard.html#airr-c-format)
#'   into an [immundata::ImmunData] object.
#'
#' Counts are taken from `duplicate_count`, or from `umi_count` if
#' `duplicate_count` is absent. If neither column is present, no count column
#' is passed to `immundata`. The selected count column must contain finite,
#' non-negative integers with no missing values.
#' @inheritParams read_airr_bulk
#' @return An [immundata::ImmunData] object containing the imported repertoires
#'   and available annotations. Columns retained by `column_profile` keep their
#'   original names, including `cell_id`, `locus`,
#'   `v_call`, `d_call`, `j_call`, `junction`, `junction_aa`,
#'   `duplicate_count`, `umi_count`, and `consensus_count`. These count columns
#'   remain separate. `umi_count` counts distinct UMIs, `duplicate_count`
#'   counts copies, and `consensus_count` counts reads supporting consensus
#'   sequences. The reader does not treat `duplicate_count` as a UMI count.
#'   Other available columns depend on the input and `column_profile`. Missing
#'   source fields remain absent. Source and manifest annotations keep their
#'   names.
#' @seealso [read_immcantation_sc()] for single-cell tables.
#' @examples
#' \dontrun{
#' path <- system.file("extdata", "immcantation", "bulk.tsv.gz",
#'   package = "immunread", mustWork = TRUE
#' )
#' idata <- read_immcantation_bulk(path,
#'   schema = c("v_call", "junction"),
#'   repertoire_schema = "sample_id", output_folder = tempfile("immunread-")
#' )
#' }
#' @export
read_immcantation_bulk <- function(path, schema, manifest = NULL,
                                   manifest_file_col = "file", column_profile = "standard",
                                   repertoire_schema = "<auto>", ...) {
  columns <- source_columns(path, manifest, manifest_file_col)
  count_col <- first_column(columns, c("duplicate_count", "umi_count"))
  prepare <- function(data, manifest = NULL) {
    profile_columns(
      keep_productive(data), schema, column_profile,
      repertoire_schema, count_col
    )
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    count_col = count_col,
    locus_col = if ("locus" %in% columns) "locus" else NULL,
    preprocess = list(profile = prepare),
    rename_columns = NULL, repertoire_schema = repertoire_schema, ...
  )
}

#' @title Read Immcantation AIRR-C single-cell tables
#' @description Reads single-cell Immcantation tables using the
#'   [AIRR-C format](https://changeo.readthedocs.io/en/stable/standard.html#airr-c-format)
#'   into an [immundata::ImmunData] object.
#'
#' The reader passes available chains to `immundata`, which selects chains
#' according to `schema`. Chains are ranked by `umi_count` when present,
#' then by `duplicate_count` if `umi_count` is absent, and finally by
#' `consensus_count` if both are absent. The selected count column must contain
#' finite, non-negative integers with no missing values. Missing values cause
#' an error even when another count column has values.
#' @inheritParams read_airr_sc
#' @inherit read_immcantation_bulk return
#' @seealso [read_immcantation_bulk()] for bulk tables.
#' @examples
#' \dontrun{
#' path <- system.file("extdata", "immcantation", "single-cell.tsv.gz",
#'   package = "immunread", mustWork = TRUE
#' )
#' idata <- read_immcantation_sc(path,
#'   schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-")
#' )
#' }
#' @export
read_immcantation_sc <- function(path, schema, manifest = NULL,
                                 manifest_file_col = "file", column_profile = "standard",
                                 repertoire_schema = "<auto>", ...) {
  columns <- source_columns(path, manifest, manifest_file_col)
  rank_col <- first_column(
    columns, c("umi_count", "duplicate_count", "consensus_count")
  )
  prepare <- function(data, manifest = NULL) {
    profile_columns(
      keep_productive(data), schema, column_profile,
      repertoire_schema, rank_col
    )
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = "cell_id", count_col = NULL,
    locus_col = if ("locus" %in% columns) "locus" else NULL,
    umi_col = rank_col, preprocess = list(profile = prepare),
    rename_columns = NULL, repertoire_schema = repertoire_schema, ...
  )
}
