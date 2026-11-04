#' @title Read Dandelion AIRR-like contigs
#' @description Reads Dandelion
#'   [AIRR contig tables](https://sc-dandelion.readthedocs.io/en/latest/notebooks/base/0_dandelion_primer.html)
#'   (the `Dandelion.data` structure) and reads their cell, locus,
#'   receptor, and abundance fields into an [immundata::ImmunData] object.
#'
#' The reader passes available chains to `immundata`, which selects chains
#' according to `schema` and count ranking. It ranks chains by `umi_count` when
#' present. If that field is absent, it uses `duplicate_count`. If both fields
#' are absent, it uses `consensus_count`. Counts used for ranking must be
#' finite whole numbers that are zero or greater, with no missing values.
#' @inheritParams read_airr_sc
#' @return An [immundata::ImmunData] object containing the imported repertoires
#'   and available annotations. **Column names kept unchanged:** Supplied AIRR
#'   names such as `cell_id`, `locus`, `v_call`, `d_call`, `j_call`, `junction`,
#'   `junction_aa`, and `consensus_count` keep their names.
#'   When `umi_count` is absent and `duplicate_count` is present, the source
#'   `duplicate_count` column is renamed to `umi_count`. Its original value is
#'   also kept as `source_duplicate_count` in the `full` column profile. When
#'   `umi_count` is already present, `duplicate_count` keeps its name and its
#'   separate copy count. `umi_count` counts distinct UMIs and
#'   `consensus_count` counts reads supporting consensus sequences. Other
#'   available columns depend on the input and `column_profile`. Source and
#'   manifest annotations keep their names.
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
  count_col <- if (is.null(rank_col)) {
    NULL
  } else if (identical(rank_col, "consensus_count")) {
    "consensus_count"
  } else {
    "umi_count"
  }
  prepare <- function(data, manifest = NULL) {
    data <- format_umi_alias(data)
    profile_columns(
      keep_productive(data), schema, column_profile, repertoire_schema,
      count_col = count_col
    )
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = "cell_id",
    locus_col = if ("locus" %in% columns) "locus" else NULL,
    umi_col = count_col, preprocess = list(profile = prepare),
    rename_columns = NULL, repertoire_schema = repertoire_schema, ...
  )
}
