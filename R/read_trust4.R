#' @title Read TRUST4 cell tables
#' @description Reads TRUST4 cell tables and separates their chain records
#'   into rows with AIRR names for the gene, sequence, and count annotations.
#'
#' The [TRUST4 barcode report format](https://github.com/liulab-dfci/TRUST4#10x-genomics-data-and-barcode-based-single-cell-data)
#' describes the embedded comma-separated chain fields. This reader expects
#' the tabular layout in the bundled example, with `X.barcode`, `chain1`,
#' `chain2`, `secondary_chain1`, and `secondary_chain2` columns.
#' See the TRUST4 paper, \doi{10.1038/s41592-021-01142-2}.
#'
#' Each non-empty primary or secondary chain record becomes a row before
#' 'immundata' selects the chains used to construct receptors.
#' @inheritParams read_airr_sc
#' @return An [immundata::ImmunData] object containing the imported repertoires
#'   and available chain annotations. The reader splits `chain1`, `chain2`,
#'   `secondary_chain1`, and `secondary_chain2` at semicolons and makes each
#'   record into its own row. Within each record, comma-separated
#'   fields 1 through 8 map in order to `v_call`, `d_call`, `j_call`, `c_call`,
#'   `junction`, `junction_aa`, `source_read_count`, and `sequence_id`.
#'   `X.barcode` maps to `cell_id`, and the field containing the record is
#'   recorded as `source_chain`. The `locus` value comes from the first three
#'   characters of `v_call`, or from `j_call` when `v_call` is `*`.
#'
#'   The read-support estimate in field 7 becomes `source_read_count`. Its
#'   fractional part is removed to give `consensus_count`, matching the
#'   upstream AIRR exporter. For example, a read-support estimate of 7.9 gives
#'   `consensus_count` 7. TRUST4 provides no separate UMI count, so
#'   `consensus_count` is used for chain ranking. Empty and `*` chain fields
#'   are skipped. With `column_profile = "full"`, original input columns such as `X.barcode` and
#'   the embedded chain fields remain under their original names, and the
#'   mapped columns are added. **These input column names are kept unchanged.**
#'   `source_read_count` is also available with this profile. The contents of
#'   the chain fields are split into the mapped columns above.
#'   With `column_profile = "minimal"` or `"standard"`, the reader
#'   keeps the columns included in that profile. The available columns also
#'   depend on the input and `schema`.
#' @examples
#' # Use one thread for this small example.
#' duckplyr::db_exec("SET threads TO 1")
#'
#' path <- system.file("extdata", "trust4", "1.1.7", "contigs.csv.gz",
#'   package = "immunread", mustWork = TRUE)
#' idata <- read_trust4(path, schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-"))
#' @export
read_trust4 <- function(path, schema, manifest = NULL,
                        manifest_file_col = "file", column_profile = "standard",
                        repertoire_schema = "<auto>", ...) {
  normalize <- function(data) {
    groups <- c("chain1", "chain2", "secondary_chain1", "secondary_chain2")
    chains <- lapply(groups, function(group) {
      data |>
        dplyr::mutate(
          cell_id = `X.barcode`, source_chain = group,
          source_record = !!rlang::sym(group)
        ) |>
        dplyr::filter(
          !is.na(source_record), source_record != "*",
          source_record != ""
        )
    })
    data <- Reduce(dplyr::union_all, chains)
    parts <- list()
    remaining <- data
    repeat {
      parts[[length(parts) + 1L]] <- dplyr::mutate(
        remaining,
        source_record = sub(";.*$", "", source_record)
      )
      more <- dplyr::filter(remaining, grepl(";", source_record))
      if (dplyr::collect(dplyr::summarise(more, n = dplyr::n()))$n[[1]] == 0L) break
      remaining <- dplyr::mutate(more,
        source_record = sub("^[^;]*;", "", source_record)
      )
    }
    data <- Reduce(dplyr::union_all, parts)
    field <- function(n) {
      if (n == 1L) {
        return(rlang::expr(sub(",.*$", "", source_record)))
      }
      pattern <- paste0("^([^,]*,){", n - 1L, "}([^,]*).*$")
      rlang::call2("sub", pattern, "\\2", rlang::sym("source_record"))
    }
    dplyr::mutate(data,
      v_call = !!field(1), d_call = !!field(2), j_call = !!field(3),
      c_call = !!field(4), junction = !!field(5),
      junction_aa = !!field(6), sequence_id = !!field(8),
      source_read_count = dd$cast_to_type(!!field(7), 0.0),
      consensus_count = dd$floor(dd$cast_to_type(!!field(7), 0.0)),
      locus = dplyr::if_else(v_call == "*", substr(j_call, 1L, 3L),
        substr(v_call, 1L, 3L)
      )
    )
  }
  prepare <- function(data, manifest = NULL) {
    profile_columns(
      keep_productive(normalize(data)), schema, column_profile,
      repertoire_schema, count_col = "consensus_count"
    )
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = "cell_id", locus_col = "locus", umi_col = "consensus_count",
    preprocess = list(trust4 = prepare), rename_columns = NULL,
    repertoire_schema = repertoire_schema, ...
  )
}
