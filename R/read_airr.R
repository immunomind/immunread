# AIRR inputs already use canonical field names. Only the productive filter and
# requested column projection are needed before immundata builds receptors.
prepare_airr <- function(data, schema, column_profile, repertoire_schema, count_col) {
  if (!"productive" %in% names(data)) {
    stop("Missing required column: productive", call. = FALSE)
  }
  data <- keep_productive(data)
  profile_columns(data, schema, column_profile, repertoire_schema, count_col)
}

#' @title Read an AIRR Rearrangement table as bulk repertoire data
#' @description Reads a bulk
#'   [AIRR Rearrangement table](https://docs.airr-community.org/en/stable/datarep/rearrangements.html).
#'
#' @param path One or more repertoire file paths or a file pattern such as
#'   `"data/*.tsv"`. Use `"<manifest>"` to read the paths from `manifest`.
#' @param schema Column names that define receptor identity, such as
#'   `c("v_call", "junction_aa")`. Use [immundata::make_receptor_schema()] to
#'   also select one chain or a pair of chains from each cell.
#' @param manifest An optional data frame with file paths and sample annotations.
#' @param manifest_file_col The name of the manifest column with file paths when
#'   `path = "<manifest>"`.
#' @param column_profile Which input annotations to keep. `"minimal"` keeps
#'   core identifiers, gene calls, junctions, and counts. `"standard"` also keeps
#'   sequence, region, and alignment annotations. `"full"` keeps all input
#'   annotations after the changes described under Value. Columns required by
#'   `schema` and `repertoire_schema` are kept with every profile.
#' @param repertoire_schema Columns that define repertoire groups. This is
#'   passed to [immundata::read_repertoires()].
#' @param ... Additional arguments passed to [immundata::read_repertoires()].
#'
#' @details Counts are taken from the first available column in this order:
#'   `umi_count`, `duplicate_count`, `consensus_count`. If none is present, no
#'   count column is passed to `immundata`.
#'
#'   So in short, this default prioritizes counts of tagged molecules over read
#'   counts. Please mind that the
#'   [AIRR schema v1.4.1 release notes](https://docs.airr-community.org/en/stable/datarep/news.html#version-1-4-1-august-27-2022)
#'   distinguish distinct UMIs (`umi_count`) from duplicate sequence observations
#'   (`duplicate_count`), and this priority order is our interpretation of the standard.
#'
#'   Count values must be finite, non-negative integers. Optional count columns
#'   can contain missing values, but the selected count column must have a value
#'   for every productive input row. Missing counts cause an error; they are not
#'   replaced with zero or one.
#'
#'   In recognized sequence and count fields, text such as `"NA"` and `"None"`
#'   becomes a missing value (`NA` in R). A row is excluded if it has no value for
#'   a receptor feature required by `schema`. These changes do not apply to
#'   other source annotations or to metadata added from `manifest`.
#' @return An [immundata::ImmunData] object containing the imported repertoire
#'   data. The chain annotations keep the supplied AIRR column names, including
#'   `v_call`, `d_call`, `j_call`, `junction`, and `junction_aa` when present.
#'
#'   `duplicate_count`, `umi_count`, and `consensus_count` remain separate
#'   columns with their original names. In AIRR terminology, they
#'   describe sequence copies, distinct UMIs, and reads used to build a consensus
#'   or contig. Selecting a count column does not rename the other count columns
#'   or change their units.
#'   The reader also keeps `cdr3` and `cdr3_aa` separate from the junction fields.
#'   It does not use one name as a replacement for another.
#'
#'   The available annotations depend on the input and `column_profile`.
#'   'immundata' adds its own analysis columns, and `manifest` can add sample
#'   annotations under their original column names.
#' @examples
#' # Use one thread for this small example.
#' duckplyr::db_exec("SET threads TO 1")
#'
#' path <- system.file("extdata", "airr", "2.0-bulk", "contigs.tsv.gz",
#'   package = "immunread", mustWork = TRUE
#' )
#' idata <- read_airr_bulk(path,
#'   schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-")
#' )
#' @export
read_airr_bulk <- function(path, schema, manifest = NULL,
                           manifest_file_col = "file", column_profile = "standard",
                           repertoire_schema = "<auto>", ...) {
  columns <- source_columns(path, manifest, manifest_file_col)
  count_col <- first_column(
    columns, c("umi_count", "duplicate_count", "consensus_count")
  )
  prepare <- function(data, manifest = NULL) {
    prepare_airr(data, schema, column_profile, repertoire_schema, count_col)
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    count_col = count_col, preprocess = list(airr = prepare),
    rename_columns = NULL, repertoire_schema = repertoire_schema, ...
  )
}

#' @title Read an AIRR Rearrangement table as single-cell repertoire data
#' @description Reads a single-cell
#'   [AIRR Rearrangement table](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)
#'   and uses its cell identifiers and counts to select receptor chains.
#'
#' To select chains, the function uses `umi_count` when it is present. If that
#' column is absent, it uses `duplicate_count`. If both are absent, it uses
#' `consensus_count`. This is the same UMI-first priority used by
#' [read_airr_bulk()], but here the counts are used to rank chains within each
#' cell. Missing values in the selected field cause an error even when another
#' count field has values. The count and missing-value rules described in
#' [read_airr_bulk()] also apply here.
#'
#' @inheritParams read_airr_bulk
#' @return An [immundata::ImmunData] object containing the imported single-cell
#'   repertoire data. If `cell_id` is absent, values from `barcode` are copied
#'   into a new `cell_id` column. If `cell_id` is already present, it is used
#'   without changing its name. When a source `barcode` column is supplied, it
#'   is also available with `column_profile = "full"`.
#'
#'   Other AIRR columns keep their original names. This includes gene calls,
#'   `locus`, `junction`, `junction_aa`, `duplicate_count`, `umi_count`, and
#'   `consensus_count` when present. The function does not rename a count field
#'   when it uses that field to select chains. These count names describe
#'   sequence copies, distinct UMIs, and reads used to build a consensus or
#'   contig in AIRR terminology. The reader also keeps `cdr3` and
#'   `cdr3_aa` separate from the junction fields.
#'
#'   The available annotations depend on the input and `column_profile`.
#'   'immundata' adds its own analysis columns, and `manifest` can add sample
#'   annotations under their original column names.
#' @examples
#' # Use one thread for this small example.
#' duckplyr::db_exec("SET threads TO 1")
#'
#' path <- system.file("extdata", "airr", "2.0-sc", "contigs.tsv.gz",
#'   package = "immunread", mustWork = TRUE
#' )
#' schema <- immundata::make_receptor_schema(
#'   c("v_call", "junction_aa"),
#'   chains = "TRB"
#' )
#' idata <- read_airr_sc(path,
#'   schema = schema,
#'   output_folder = tempfile("immunread-")
#' )
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
    prepare_airr(data, schema, column_profile, repertoire_schema, rank_col)
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = "cell_id", locus_col = if ("locus" %in% columns) "locus" else NULL,
    umi_col = rank_col, preprocess = list(airr = prepare),
    rename_columns = NULL, repertoire_schema = repertoire_schema, ...
  )
}
