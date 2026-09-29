test_that("AIRR 2.0-sc has expected chains and receptors", {
  idata <- read_fixture(read_airr_sc, c("airr", "2.0-sc", "contigs.tsv.gz"))
  expect_s3_class(idata, "ImmunData")
  expect_true(all(c("cell_id", "locus", "v_call", "junction_aa") %in% names(annotations(idata))))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 1L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 1L)
})

test_that("AIRR single-cell ranks same-locus chains by umi_count", {
  idata <- read_fixture(read_airr_sc, c("airr", "2.0-sc", "contigs.tsv.gz"))
  a <- annotations(idata)
  selected <- a$junction_aa[a$cell_id == "cell-a"]
  expect_equal(selected, "CASRF")
  expect_false("cell-b" %in% a$cell_id)
})

test_that("AIRR single-cell falls back to duplicate then consensus counts", {
  read_ranked <- function(extra) {
    path <- tempfile(fileext = ".tsv")
    x <- data.frame(
      sequence_id = c("rank-1", "rank-2"), junction = c("TGTGCT", "TGTAGT"),
      d_call = "", cell_id = c("cell-a", "cell-a"), locus = "TRB", productive = "T",
      v_call = "TRBV1*01", junction_aa = c("CASSF", "CASRF"),
      consensus_count = extra$consensus_count
    )
    if (!is.null(extra$duplicate_count)) x$duplicate_count <- extra$duplicate_count
    utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
    idata <- read_airr_sc(path, c("v_call", "junction_aa"), output_folder = tempfile())
    a <- annotations(idata)
    a$junction_aa[a$cell_id == "cell-a"]
  }

  expect_equal(read_ranked(list(duplicate_count = c(2, 8), consensus_count = c(9, 1))), "CASRF")
  expect_equal(read_ranked(list(duplicate_count = NULL, consensus_count = c(2, 8))), "CASRF")
})

test_that("bulk-derived 1.6-sc fixture cannot substitute barcodes for ranking counts", {
  path <- fixture_path("airr", "1.6-sc", "contigs.tsv.gz")
  expect_error(read_airr_sc(path, c("v_call", "cdr3_aa"), output_folder = tempfile()),
               "umi_col")
})

test_that("AIRR single-cell accepts barcode and prefers an existing cell_id", {
  path <- tempfile(fileext = ".tsv")
  data <- data.frame(barcode = "cell-a", productive = "T", locus = "TRB",
                     v_call = "TRBV1*01", junction_aa = "CASSF", umi_count = 2L)
  for (has_cell_id in c(FALSE, TRUE)) {
    if (has_cell_id) data$cell_id <- "cell-b"
    utils::write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE)
    idata <- read_airr_sc(path, c("v_call", "junction_aa"), output_folder = tempfile())
    expect_equal(annotations(idata)$cell_id, if (has_cell_id) "cell-b" else "cell-a")
    expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 1L)
  }
})

test_that("AIRR single-cell requires a real UMI or read count", {
  path <- tempfile(fileext = ".tsv")
  utils::write.table(data.frame(
    cell_id = "cell-a", locus = "TRB", productive = "T",
    v_call = "TRBV1*01", junction_aa = "CASSF"
  ), path, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(read_airr_sc(path, c("v_call", "junction_aa"),
    output_folder = tempfile()), "requires.*umi_col")
})
