expected <- list(`1.6-bulk` = c(chains = 80L, receptors = 67L),
                 `2.0-bulk` = c(chains = 2L, receptors = 2L))
for (version in names(expected)) {
  test_that(paste("AIRR", version, "has expected chains and receptors"), {
    idata <- read_fixture(read_airr_bulk, c("airr", version, "contigs.tsv.gz"))
    expect_s3_class(idata, "ImmunData")
    expect_true(all(c("v_call", "junction_aa") %in% names(annotations(idata))))
    expect_equal(idata |> dplyr::count() |> dplyr::pull(n), expected[[version]][["chains"]])
    expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), expected[[version]][["receptors"]])
    if (version == "2.0-bulk") {
      expect_true(all(c("duplicate_count", "umi_count") %in% names(annotations(idata))))
    }
  })
}

test_that("bulk-derived 1.6-sc fixture has expected chains and receptors in bulk mode", {
  # This source has no single-cell ranking counts; its barcodes are synthetic.
  # cdr3_aa is populated, whereas junction_aa is empty.
  idata <- read_fixture(read_airr_bulk, c("airr", "1.6-sc", "contigs.tsv.gz"),
                        schema = c("v_call", "cdr3_aa"))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 802L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 737L)
})

test_that("AIRR bulk uses duplicate_count before umi_count and keeps count meanings", {
  idata <- read_fixture(read_airr_bulk, c("airr", "2.0-bulk", "contigs.tsv.gz"))
  a <- annotations(idata)
  expect_equal(a$imd_n_chains, a$duplicate_count)
  expect_false(identical(a$duplicate_count, a$umi_count))
})

test_that("AIRR column profiles prune source columns while protecting analysis fields", {
  path <- fixture_path("airr", "2.0-bulk", "contigs.tsv.gz")
  minimal <- read_airr_bulk(path, c("v_call", "junction_aa"), column_profile = "minimal", output_folder = tempfile())
  standard <- read_airr_bulk(path, c("v_call", "junction_aa"), column_profile = "standard", output_folder = tempfile())
  full <- read_airr_bulk(path, c("v_call", "junction_aa"), column_profile = "full", output_folder = tempfile())
  min_names <- names(annotations(minimal))
  std_names <- names(annotations(standard))
  full_names <- names(annotations(full))

  expect_true(all(c("v_call", "junction_aa", "productive", "duplicate_count", "repertoire_id", "sample_id", "subject_id") %in% min_names))
  expect_false("custom_note" %in% min_names)
  expect_true("sequence" %in% std_names)
  expect_false("custom_note" %in% std_names)
  expect_true(all(c("custom_note", "sequence_alignment", "umi_count", "consensus_count") %in% full_names))

  custom_schema <- read_airr_bulk(path, c("v_call", "junction_aa", "custom_note"),
    column_profile = "minimal", output_folder = tempfile())
  expect_true("custom_note" %in% names(annotations(custom_schema)))

  source <- tempfile(fileext = ".tsv")
  utils::write.table(data.frame(
    sequence_id = "seq-1", productive = "T", v_call = "TRBV1*01",
    junction_aa = "CASSF", duplicate_count = 1L,
    sequence_alignment = "ACGT", germline_alignment = "ACGT",
    v_sequence_start = 1L, custom_alignment = "ACGT", custom_start = 1L
  ), source, sep = "\t", quote = FALSE, row.names = FALSE)
  standard_custom <- read_airr_bulk(source, c("v_call", "junction_aa"),
    output_folder = tempfile()) |>
    annotations()
  expect_true(all(c("sequence_alignment", "germline_alignment", "v_sequence_start") %in% names(standard_custom)))
  expect_false(any(c("custom_alignment", "custom_start") %in% names(standard_custom)))
})

test_that("AIRR validation rejects missing required columns", {
  missing_productive <- tempfile(fileext = ".tsv")
  writeLines(c("v_call\tjunction_aa", "TRBV1*01\tCASSF"), missing_productive)
  expect_error(read_airr_bulk(missing_productive, c("v_call", "junction_aa"), output_folder = tempfile()), "productive")

  missing_feature <- tempfile(fileext = ".tsv")
  writeLines(c("productive\tv_call", "T\tTRBV1*01"), missing_feature)
  expect_error(read_airr_bulk(missing_feature, c("v_call", "junction_aa"), output_folder = tempfile()), "junction_aa")
})

test_that("manifest columns annotate AIRR input", {
  path <- fixture_path("airr", "2.0-bulk", "contigs.tsv.gz")
  manifest <- data.frame(file = path, study_label = "keep me")
  idata <- read_airr_bulk("<manifest>", c("v_call", "junction_aa"),
    manifest = manifest, manifest_file_col = "file", column_profile = "minimal",
    output_folder = tempfile())

  expect_true("study_label" %in% names(annotations(idata)))
  expect_true(all(annotations(idata)$study_label == "keep me"))
})
