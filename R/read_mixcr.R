#' @title Read MiXCR single-cell clones
#' @description Reads MiXCR clone tables and maps their cell, chain, gene,
#'   sequence, and count fields to AIRR-style columns.
#'
#' Supplied secondary chains are retained; `immundata` ranks chains when
#' constructing receptors. Molecule and read counts remain distinct.
#' @inheritParams read_airr_sc
#' @examples
#' path <- system.file("extdata", "mixcr", "4.7.0", "contigs.csv.gz",
#'   package = "immunread", mustWork = TRUE)
#' idata <- read_mixcr(path, schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-"))
#' @export
read_mixcr <- function(path, schema, manifest = NULL,
                       manifest_file_col = "file", column_profile = "standard",
                       repertoire_schema = "<auto>", ...) {
  prepare <- function(data, manifest = NULL) {
    data <- format_aliases(data, c(
      cell_id = "cellGroup", locus = "topChains",
      v_call = "allVHitsWithScore", d_call = "allDHitsWithScore",
      j_call = "allJHitsWithScore", junction = "nSeqCDR3",
      junction_aa = "aaSeqCDR3", umi_count = "uniqueMoleculeCount",
      duplicate_count = "readCount"
    ))
    for (name in intersect(c("v_call", "d_call", "j_call"), names(data))) {
      data <- dplyr::mutate(
        data,
        !!name := gsub("\\([^)]*\\)", "", !!rlang::sym(name))
      )
    }
    profile_columns(
      keep_productive(data), schema, column_profile,
      repertoire_schema
    )
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = "cell_id", locus_col = "locus", umi_col = "umi_count",
    preprocess = list(mixcr = prepare), rename_columns = NULL,
    repertoire_schema = repertoire_schema, ...
  )
}
