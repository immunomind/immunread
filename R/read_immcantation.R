#' @title Read Immcantation AIRR-like tables
#' @description Reads Immcantation tables using the
#'   [AIRR-C format](https://changeo.readthedocs.io/en/stable/standard.html#airr-c-format)
#'   in bulk or single-cell mode and reads their annotations into an
#'   [immundata::ImmunData] object.
#'
#' In single-cell mode, the reader passes available chains to `immundata`,
#' which selects chains according to `schema` and count ranking. In bulk mode,
#' abundance uses `duplicate_count` when that field is present. If it is
#' absent, abundance uses `umi_count`. In single-cell mode, chains are ranked
#' by `umi_count` when present. If it is absent, ranking uses
#' `duplicate_count`. If both fields are absent, ranking uses
#' `consensus_count`. Counts used for abundance or ranking must be finite whole
#' numbers that are zero or greater, with no missing values. Unavailable
#' source annotations remain missing and are not reconstructed from other
#' regions.
#' @inheritParams read_airr_sc
#' @return An [immundata::ImmunData] object containing the imported repertoires
#'   and available annotations. **Column names kept unchanged:** The reader
#'   does not convert input column names. AIRR fields kept by the selected
#'   `column_profile` keep their supplied names, including `cell_id`, `locus`,
#'   `v_call`, `d_call`, `j_call`, `junction`, `junction_aa`,
#'   `duplicate_count`, `umi_count`, and `consensus_count`. These count columns
#'   remain separate. `umi_count` counts distinct UMIs, `duplicate_count`
#'   counts copies, and `consensus_count` counts reads supporting consensus
#'   sequences. The reader does not treat `duplicate_count` as a UMI count.
#'   Other available columns depend on the input and `column_profile`. Missing
#'   source fields remain absent. Source and manifest annotations keep their
#'   names.
#' @param mode `"sc"` for single-cell tables or `"bulk"` for bulk tables.
#' @examples
#' path <- system.file("extdata", "immcantation", "single-cell.tsv.gz",
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
                    repertoire_schema, if (mode == "bulk") count_col else rank_col)
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
