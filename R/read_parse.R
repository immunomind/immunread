#' @title Read Parse Biosciences TCR tables
#' @description Reads wide Parse Biosciences TCR tables, emits one row per
#'   supplied chain, and maps the result to AIRR-style columns.
#'
#' Each supplied primary or secondary TRA/TRB chain becomes a row. `immundata`
#' ranks chains when constructing receptors. Transcript counts are used for
#' ranking; source read counts remain separate.
#' @inheritParams read_airr_sc
#' @examples
#' path <- system.file("extdata", "parsebio", "1.6.1", "contigs.csv.gz",
#'   package = "immunread", mustWork = TRUE)
#' idata <- read_parsebio(path, schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-"))
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
    chains <- lapply(groups, function(group) {
      aliases <- c(
        cell_id = "Barcode", v_call = paste0(group, "_V"),
        d_call = paste0(group, "_D"),
        j_call = paste0(group, "_J"),
        c_call = paste0(group, "_C"),
        junction_aa = paste0(group, "_cdr3_aa"),
        umi_count = paste0(group, "_transcript_count")
      )
      format_chain_aliases(data, aliases, targets) |>
        dplyr::mutate(umi_count = as.integer(umi_count)) |>
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
      repertoire_schema
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
