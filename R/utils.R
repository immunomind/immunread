# Fields retained by the compact input profiles. User manifest columns are joined
# downstream, after this projection.
.minimal_columns <- c(
  "sequence_id", "repertoire_id", "sample_id", "subject_id",
  "sample_processing_id", "data_processing_id", "cell_id", "clone_id",
  "productive", "locus", "locus_species", "v_call", "d_call", "d2_call",
  "j_call", "c_call", "junction", "junction_aa", "junction_length",
  "junction_aa_length", "duplicate_count", "umi_count", "consensus_count"
)

.standard_columns <- c(
  .minimal_columns, "sequence", "sequence_aa", "rev_comp", "vj_in_frame",
  "stop_codon", "complete_vdj", "v_frameshift", "j_frameshift",
  "d_frame", "d2_frame", "rearrangement_type", "reactivity_id",
  "reactivity_ref", paste0(rep(
    c(
      "fwr1", "cdr1", "fwr2", "cdr2",
      "fwr3", "cdr3", "fwr4", "np1", "np2", "np3"
    ),
    each = 2
  ), c("", "_aa")),
  "np1_length", "np2_length", "np3_length", "n1_length", "n2_length",
  "n3_length", "p3v_length", "p5d_length", "p3d_length", "p5d2_length",
  "p3d2_length", "p5j_length"
)

# immundata needs the abundance/rank column name before preprocessing runs.
# This inspects the first input schema lazily; immundata validates the files.
source_columns <- function(path, manifest, manifest_file_col) {
  file <- if (identical(path, "<manifest>")) {
    manifest[[manifest_file_col]][1]
  } else {
    Sys.glob(path)[1]
  }
  type <- tolower(tools::file_ext(sub("\\.(gz|gzip)$", "", file)))
  if (type == "parquet") {
    names(duckplyr::read_parquet_duckdb(file, prudence = "stingy"))
  } else {
    names(duckplyr::read_csv_duckdb(
      file,
      prudence = "stingy",
      options = list(delim = if (type == "tsv") "\t" else ",")
    ))
  }
}

first_column <- function(columns, choices) {
  found <- intersect(choices, columns)
  if (length(found)) found[1] else NULL
}

keep_productive <- function(data) {
  if (!"productive" %in% names(data)) {
    return(data)
  }
  immundata::make_productive_filter(
    truthy = c(
      "TRUE", "True", "true", "T", "t", "YES", "Yes", "yes",
      "Y", "y", "1"
    )
  )(data)
}

profile_columns <- function(data, schema, column_profile, repertoire_schema) {
  column_profile <- match.arg(column_profile, c("minimal", "standard", "full"))
  if (column_profile == "full" || is.function(repertoire_schema)) {
    return(data)
  }

  features <- if (is.character(schema)) schema else schema$features
  protected <- c(
    features, "imd_filename", "productive", "cell_id", "locus",
    "duplicate_count", "umi_count", "consensus_count"
  )
  if (is.character(repertoire_schema) &&
    !identical(repertoire_schema, "<auto>") &&
    !identical(repertoire_schema, "<manifest>")) {
    protected <- c(protected, repertoire_schema)
  }
  keep <- if (column_profile == "minimal") {
    .minimal_columns
  } else {
    c(
      .standard_columns,
      grep(
        "^(sequence_alignment|germline_alignment|(v|d|d2|j|c)_(alignment|sequence_alignment|germline_alignment)|(v|d|d2|j|c)_(sequence|germline)_(start|end)|(fwr[1-4]|cdr[1-3]|junction)_(start|end))$",
        names(data),
        value = TRUE
      )
    )
  }
  dplyr::select(data, dplyr::any_of(unique(c(keep, protected))))
}

# Add canonical fields while keeping source columns available to the full profile.
format_aliases <- function(data, aliases) {
  aliases <- aliases[unname(aliases) %in% names(data) &
    !names(aliases) %in% names(data)]
  if (!length(aliases)) {
    return(data)
  }
  dplyr::mutate(data, !!!stats::setNames(rlang::syms(unname(aliases)), names(aliases)))
}

# Align emitted chain rows without inventing an abundance when none was supplied.
format_chain_aliases <- function(data, aliases, targets) {
  data <- format_aliases(data, aliases)
  missing <- setdiff(targets, names(data))
  if ("umi_count" %in% missing) {
    stop("Missing per-chain UMI/read count: ", aliases[["umi_count"]],
         call. = FALSE)
  }
  if (!length(missing)) {
    return(data)
  }
  values <- rep(list(NA_character_), length(missing))
  dplyr::mutate(data, !!!stats::setNames(values, missing))
}
