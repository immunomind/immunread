#' @title Read WAT3R TCR tables
#' @description Reads wide WAT3R TCR tables, emits one row per supplied chain,
#'   and maps their receptor and UMI fields to AIRR-style columns.
#'
#' Each supplied TRA, TRB, or secondary TRA chain becomes a row. `immundata`
#' ranks chains when constructing receptors. UMI and read counts remain separate.
#' @inheritParams read_airr_sc
#' @examples
#' path <- system.file("extdata", "wat3r", "1.3", "contigs.csv.gz",
#'   package = "immunread", mustWork = TRUE)
#' idata <- read_wat3r(path, schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-"))
#' @export
read_wat3r <- function(path, schema, manifest = NULL,
                       manifest_file_col = "file", column_profile = "standard",
                       repertoire_schema = "<auto>", ...) {
  normalize <- function(data) {
    # duckplyr union_all cannot combine the source's logical InRNAseq column.
    if ("InRNAseq" %in% names(data)) {
      data <- dplyr::mutate(data, InRNAseq = dplyr::if_else(
        InRNAseq, "TRUE", "FALSE"))
    }
    groups <- c("TRA", "TRB", "TRA.2")
    groups <- groups[paste0(groups, "_CDR3") %in% names(data)]
    targets <- c("cell_id", "v_call", "d_call", "j_call",
                 "junction", "junction_aa")
    if (any(paste0(groups, "_CDR3_UMIcount") %in% names(data))) {
      targets <- c(targets, "umi_count")
    }
    chains <- lapply(groups, function(group) {
      suffix <- if (group == "TRA.2") ".2" else ""
      locus <- if (group == "TRA.2") "TRA" else group
      aliases <- c(cell_id = "BC",
                   v_call = paste0(locus, "V", suffix),
                   j_call = paste0(locus, "J", suffix),
                   junction = paste0(group, "_CDR3nuc"),
                   junction_aa = paste0(group, "_CDR3"),
                   umi_count = paste0(group, "_CDR3_UMIcount"))
      if (group == "TRB") aliases <- c(aliases, d_call = "TRBD")
      format_chain_aliases(data, aliases, targets) |>
        dplyr::mutate(umi_count = as.integer(umi_count)) |>
        dplyr::mutate(locus = locus, source_chain = group) |>
        dplyr::filter(!is.na(junction_aa), junction_aa != "", junction_aa != "NA")
    })
    Reduce(dplyr::union_all, chains)
  }
  prepare <- function(data, manifest = NULL) {
    profile_columns(keep_productive(normalize(data)), schema, column_profile,
                    repertoire_schema)
  }
  immundata::read_repertoires(
    path, schema, manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = "cell_id", locus_col = "locus", umi_col = "umi_count",
    preprocess = list(wat3r = prepare), rename_columns = NULL,
    repertoire_schema = repertoire_schema, ...
  )
}
