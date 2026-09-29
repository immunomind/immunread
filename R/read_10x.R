#' @title Read 10x Genomics V(D)J contigs
#' @description Reads 10x Genomics V(D)J contig annotations and maps their
#'   barcode, locus, gene, sequence, and count fields to AIRR-style columns.
#'
#' Supplied secondary chains are retained; `immundata` ranks chains when
#' constructing receptors. `umis` and `reads` remain distinct counts.
#' @inheritParams read_airr_sc
#' @examples
#' path <- system.file("extdata", "10x", "9.0.1", "contigs.csv.gz",
#'   package = "immunread", mustWork = TRUE)
#' idata <- read_10x(path, schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-"))
#' @export
read_10x <- function(path, schema, manifest = NULL,
                     manifest_file_col = "file", column_profile = "standard",
                     repertoire_schema = "<auto>", ...) {
  prepare <- function(data, manifest = NULL) {
    data <- format_aliases(data, c(
      cell_id = "barcode", sequence_id = "contig_id", locus = "chain",
      v_call = "v_gene", d_call = "d_gene", j_call = "j_gene",
      c_call = "c_gene", junction = "cdr3_nt",
      junction_aa = "cdr3", umi_count = "umis", consensus_count = "reads"
    ))
    # 10x's unsuffixed region fields are amino acids; AIRR's are nucleotides.
    # Very incovenient, but what can we do...
    for (region in c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "fwr4")) {
      if (!region %in% names(data)) next
      data <- dplyr::mutate(data, !!paste0(region, "_aa") := !!rlang::sym(region))
      nt <- paste0(region, "_nt")
      if (nt %in% names(data)) {
        data <- dplyr::mutate(data, !!region := !!rlang::sym(nt))
      } else {
        data <- dplyr::select(data, -dplyr::all_of(region))
      }
    }
    if ("cdr3" %in% names(data)) data <- dplyr::rename(data, source_cdr3 = cdr3)
    if ("cdr3_nt" %in% names(data)) {
      data <- dplyr::rename(data, source_cdr3_nt = cdr3_nt)
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
    preprocess = list(`10x` = prepare), rename_columns = NULL,
    repertoire_schema = repertoire_schema, ...
  )
}
