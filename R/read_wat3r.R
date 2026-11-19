#' @title Read WAT3R TCR tables
#' @description Reads WAT3R
#'   [barcode_results.csv tables](https://github.com/mainciburu/WAT3R#downstream)
#'   in which each row describes one cell. The reader separates the chains into
#'   rows and uses AIRR names for the output annotations.
#'
#' The input contains `BC`, `TRA_*`, `TRB_*`, and optional `TRA.2_*` columns.
#' See the WAT3R paper, \doi{10.1093/bioinformatics/btac382}.
#'
#' The reader prepares one row for each supplied TRA, TRB, or secondary TRA
#' chain. Rows without a CDR3 amino acid sequence are excluded before
#' 'immundata' selects the chains used to construct receptors.
#' @inheritParams read_airr_sc
#' @return An [immundata::ImmunData] object containing the imported repertoires
#'   and available chain annotations. Gene column names are converted for each
#'   chain group as follows:
#'
#'   * For TRA, `TRAV` becomes `v_call` and `TRAJ` becomes `j_call`.
#'   * For TRB, `TRBV` becomes `v_call`, `TRBD` becomes `d_call`, and `TRBJ`
#'     becomes `j_call`.
#'   * For secondary TRA, `TRAV.2` becomes `v_call` and `TRAJ.2` becomes `j_call`.
#'
#'   The matching `TRA_CDR3nuc`, `TRB_CDR3nuc`, and `TRA.2_CDR3nuc` values map
#'   to `cdr3`. These nucleotide sequences exclude conserved anchor codons.
#'   The matching `*_CDR3` amino acid sequences map to `junction_aa` and include
#'   the conserved amino acids. WAT3R does not supply a nucleotide `junction`,
#'   so this reader does not create one. The matching `*_nReads` values map to
#'   `consensus_count`. `*_CDR3_UMIcount` values map to `umi_count` for
#'   chain ranking. `BC` maps to `cell_id`. The chain group supplies `locus`
#'   (`TRA.2` is recorded as `TRA`) and `source_chain`.
#'
#'   AIRR output columns that are already present keep their names and are used
#'   instead of the corresponding source columns. The reader always sets
#'   `locus` and `source_chain` from the chain group.
#'
#'   If no group has an `nReads` column,
#'   `consensus_count` is absent. If some groups have that column, chains from
#'   groups without it have `NA` for `consensus_count`. Missing optional gene
#'   calls also have `NA`.
#'
#'   With `column_profile = "full"`, the reader retains the original columns,
#'   including `BC` and the wide chain columns, and adds the AIRR columns.
#'   With `column_profile = "minimal"` or `"standard"`, the reader keeps
#'   the columns included in that profile. The available columns also depend
#'   on the input and `schema`.
#' @examples
#' # Use one thread for this small example.
#' duckplyr::db_exec("SET threads TO 1")
#'
#' path <- system.file("extdata", "wat3r", "1.3", "contigs.csv.gz",
#'   package = "immunread", mustWork = TRUE
#' )
#' idata <- read_wat3r(path,
#'   schema = c("v_call", "junction_aa"),
#'   output_folder = tempfile("immunread-")
#' )
#' @export
read_wat3r <- function(path, schema, manifest = NULL,
                       manifest_file_col = "file", column_profile = "standard",
                       repertoire_schema = "<auto>", ...) {
  normalize <- function(data) {
    # duckplyr union_all cannot combine the source's logical InRNAseq column.
    if ("InRNAseq" %in% names(data)) {
      data <- dplyr::mutate(data, InRNAseq = dplyr::if_else(
        InRNAseq, "TRUE", "FALSE"
      ))
    }
    groups <- c("TRA", "TRB", "TRA.2")
    groups <- groups[paste0(groups, "_CDR3") %in% names(data)]
    targets <- c(
      "cell_id", "v_call", "d_call", "j_call",
      "cdr3", "junction_aa"
    )
    if (any(paste0(groups, "_nReads") %in% names(data))) {
      targets <- c(targets, "consensus_count")
    }
    if (any(paste0(groups, "_CDR3_UMIcount") %in% names(data))) {
      targets <- c(targets, "umi_count")
    }
    chains <- lapply(groups, function(group) {
      suffix <- if (group == "TRA.2") ".2" else ""
      locus <- if (group == "TRA.2") "TRA" else group
      aliases <- c(
        cell_id = "BC",
        v_call = paste0(locus, "V", suffix),
        j_call = paste0(locus, "J", suffix),
        cdr3 = paste0(group, "_CDR3nuc"),
        junction_aa = paste0(group, "_CDR3"),
        umi_count = paste0(group, "_CDR3_UMIcount"),
        consensus_count = paste0(group, "_nReads")
      )
      if (group == "TRB") aliases <- c(aliases, d_call = "TRBD")
      format_chain_aliases(data, aliases, targets) |>
        dplyr::mutate(locus = locus, source_chain = group) |>
        dplyr::filter(!is.na(junction_aa), junction_aa != "", junction_aa != "NA")
    })
    Reduce(dplyr::union_all, chains)
  }
  prepare <- function(data, manifest = NULL) {
    profile_columns(keep_productive(normalize(data)), schema, column_profile,
      repertoire_schema,
      count_col = "umi_count"
    )
  }
  immundata::read_repertoires(
    path, schema,
    manifest = manifest, manifest_file_col = manifest_file_col,
    barcode_col = "cell_id", locus_col = "locus", umi_col = "umi_count",
    preprocess = list(wat3r = prepare), rename_columns = NULL,
    repertoire_schema = repertoire_schema, ...
  )
}
