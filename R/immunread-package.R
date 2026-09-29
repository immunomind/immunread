#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom immundata read_repertoires
## usethis namespace: end
NULL

# Column names used in lazy duckplyr expressions.
utils::globalVariables(c(
  ":=", "cdr3", "cdr3_nt", "frame_type", "InRNAseq", "j_call",
  "junction_aa", "sequenceStatus", "source_record", "umi_count", "v_call", "X.barcode"
))
