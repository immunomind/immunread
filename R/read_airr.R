# AIRR inputs already use canonical field names. Only the productive filter and
# requested column projection are needed before immundata builds receptors.
prepare_airr <- function(data, schema, column_profile, repertoire_schema) {
  if (!"productive" %in% names(data)) {
    stop("Missing required column: productive", call. = FALSE)
  }
  data <- keep_productive(data)
  profile_columns(data, schema, column_profile, repertoire_schema)
}

#' @title Read an AIRR Rearrangement table as bulk repertoire data
#' @description Reads a bulk AIRR Rearrangement table, filters productive
#'   rows, and passes its canonical AIRR fields to `immundata`.
#'
#' @param path Repertoire file paths, glob, or `"<manifest>"`.
#' @param schema Receptor features or an [immundata::make_receptor_schema()]
#'   result.
#' @param manifest Optional native immundata manifest.
#' @param manifest_file_col Manifest column containing file paths when
#'   `path = "<manifest>"`.
#' @param column_profile Source columns to keep: `"minimal"`, `"standard"`,
#'   or `"full"`.
#' @param repertoire_schema Repertoire grouping passed to
#'   [immundata::read_repertoires()].
#' @param ... Other arguments passed to [immundata::read_repertoires()].
#' @return An `ImmunData` object.
#' @examples
#' path <- system.file("extdata", "airr", "2.0-bulk", "contigs.tsv.gz",
#'   package = "immunread", mustWork = TRUE)
#' idata <- read_airr_bulk(path, schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-"))
#' @export
read_airr_bulk <- function(path, schema, manifest = NULL,
                           manifest_file_col = "file", column_profile = "standard",
                           repertoire_schema = "<auto>", ...) {
  columns <- source_columns(path, manifest, manifest_file_col)
  count_col <- first_column(columns, c("duplicate_count", "umi_count"))
  prepare <- function(data, manifest = NULL) {
    prepare_airr(data, schema, column_profile, repertoire_schema)
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    count_col = count_col, preprocess = list(airr = prepare),
    rename_columns = NULL, repertoire_schema = repertoire_schema, ...
  )
}

#' @title Read an AIRR Rearrangement table as single-cell repertoire data
#' @description Reads a single-cell AIRR Rearrangement table, maps `barcode`
#'   to `cell_id` when needed, and selects chains using their count field.
#'
#' Uses `barcode` as `cell_id` when `cell_id` is absent.
#'
#' @inheritParams read_airr_bulk
#' @return An `ImmunData` object.
#' @examples
#' path <- system.file("extdata", "airr", "2.0-sc", "contigs.tsv.gz",
#'   package = "immunread", mustWork = TRUE)
#' schema <- immundata::make_receptor_schema(
#'   c("v_call", "junction_aa"), chains = "TRB")
#' idata <- read_airr_sc(path, schema = schema,
#'   output_folder = tempfile("immunread-"))
#' @export
read_airr_sc <- function(path, schema, manifest = NULL,
                         manifest_file_col = "file", column_profile = "standard",
                         repertoire_schema = "<auto>", ...) {
  columns <- source_columns(path, manifest, manifest_file_col)
  rank_col <- first_column(
    columns,
    c("umi_count", "duplicate_count", "consensus_count")
  )
  prepare <- function(data, manifest = NULL) {
    data <- format_aliases(data, c(cell_id = "barcode"))
    prepare_airr(data, schema, column_profile, repertoire_schema)
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = "cell_id", locus_col = if ("locus" %in% columns) "locus" else NULL,
    umi_col = rank_col, preprocess = list(airr = prepare),
    rename_columns = NULL, repertoire_schema = repertoire_schema, ...
  )
}
