#' @title Read TRUST4 cell tables
#' @description Reads TRUST4 cell tables, expands embedded chain records, and
#'   maps their receptor and count fields to AIRR-style columns.
#'
#' Each supplied primary or secondary chain becomes a row. `immundata` ranks
#' chains when constructing receptors. The embedded `read_cnt` becomes
#' `consensus_count`, following TRUST4's AIRR export.
#' @inheritParams read_airr_sc
#' @examples
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
      consensus_count = as.integer(!!field(7)),
      locus = dplyr::if_else(v_call == "*", substr(j_call, 1L, 3L),
        substr(v_call, 1L, 3L)
      )
    )
  }
  prepare <- function(data, manifest = NULL) {
    profile_columns(
      keep_productive(normalize(data)), schema, column_profile,
      repertoire_schema
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
