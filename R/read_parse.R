#' @title Read Parse Biosciences TCR tables
#' @description Reads Parse Biosciences TCR tables in which each row describes
#'   one cell. The reader separates the chains into rows and uses AIRR names
#'   for the output annotations.
#'
#' See Parse's
#' [public immune annotation overview](https://resources.parsebiosciences.com/support-suite/hc/en-us/articles/44260645430676-understanding-the-anndata-scanpy-object-that-is-downloadable-from-trailmaker-insights-module)
#' for descriptions of the wide per-cell fields. This reader
#' expects the wide layout with `Barcode`, `TRA_*`, `TRB_*`, and optional
#' `secondary_TRA_*` / `secondary_TRB_*` columns. For AIRR-formatted
#' `tcr_annotation_airr.tsv` exports, use [read_airr_sc()].
#'
#' The reader prepares one row for each supplied primary or secondary TRA or
#' TRB chain. Rows without a CDR3 amino acid sequence are excluded before
#' 'immundata' selects the chains used to construct receptors.
#' @inheritParams read_airr_sc
#' @return An [immundata::ImmunData] object containing the imported repertoires
#'   and available chain annotations. For each of `TRA`, `TRB`,
#'   `secondary_TRA`, and `secondary_TRB`, the reader copies `*_V` to `v_call`,
#'   `*_D` to `d_call`, `*_J` to `j_call`, and `*_C` to `c_call`.
#'   It copies `*_cdr3_aa` to `junction_aa` and `*_transcript_count` to `umi_count` for chain
#'   ranking and `*_read_count` to `consensus_count`. `Barcode` maps to
#'   `cell_id`, and the group name supplies `locus` after removing the
#'   `secondary_` prefix. The group name is also recorded in `source_chain`.
#'
#'   AIRR output columns that are already present keep their names and are used
#'   instead of the corresponding source columns. The reader always sets
#'   `locus` and `source_chain` from the chain group.
#'
#'   If no group has a read-count column,
#'   `consensus_count` is absent. If some groups have that column, chains from
#'   groups without it have `NA` for `consensus_count`. Missing optional gene
#'   calls also have `NA`.
#'
#'   With `column_profile = "full"`, original source columns,
#'   including `Barcode` and the wide chain columns, remain under their original
#'   names and the mapped columns are added. **These source column names are
#'   kept unchanged.** With `column_profile = "minimal"` or `"standard"`, the
#'   reader keeps the columns included in that profile. The available columns
#'   also depend on the input and `schema`.
#' @examples
#' # Use one thread for this small example.
#' duckplyr::db_exec("SET threads TO 1")
#'
#' path <- system.file("extdata", "parsebio", "1.6.1", "contigs.csv.gz",
#'   package = "immunread", mustWork = TRUE
#' )
#' idata <- read_parsebio(path,
#'   schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-")
#' )
#' @export
read_parsebio <- function(path, schema, manifest = NULL,
                          manifest_file_col = "file", column_profile = "standard",
                          repertoire_schema = "<auto>", ...) {
  normalize <- function(data) {
    groups <- c("TRA", "TRB", "secondary_TRA", "secondary_TRB")
    groups <- groups[paste0(groups, "_cdr3_aa") %in% names(data)]
    targets <- c("cell_id", "v_call", "d_call", "j_call", "c_call", "junction_aa")
    if (any(paste0(groups, "_transcript_count") %in% names(data))) {
      targets <- c(targets, "umi_count")
    }
    if (any(paste0(groups, "_read_count") %in% names(data))) {
      targets <- c(targets, "consensus_count")
    }
    chains <- lapply(groups, function(group) {
      aliases <- c(
        cell_id = "Barcode", v_call = paste0(group, "_V"),
        d_call = paste0(group, "_D"),
        j_call = paste0(group, "_J"),
        c_call = paste0(group, "_C"),
        junction_aa = paste0(group, "_cdr3_aa"),
        umi_count = paste0(group, "_transcript_count"),
        consensus_count = paste0(group, "_read_count")
      )
      format_chain_aliases(data, aliases, targets) |>
        dplyr::mutate(
          locus = sub("^secondary_", "", group),
          source_chain = group
        ) |>
        dplyr::filter(!is.na(junction_aa), junction_aa != "", junction_aa != "NA")
    })
    Reduce(dplyr::union_all, chains)
  }
  prepare <- function(data, manifest = NULL) {
    profile_columns(
      keep_productive(normalize(data)), schema, column_profile,
      repertoire_schema,
      count_col = "umi_count"
    )
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = "cell_id", locus_col = "locus", umi_col = "umi_count",
    preprocess = list(parsebio = prepare), rename_columns = NULL,
    repertoire_schema = repertoire_schema, ...
  )
}
