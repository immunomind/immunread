test_that("reserved missing markers do not become receptor features or metadata edits", {
  path <- tempfile(fileext = ".tsv")
  utils::write.table(data.frame(
    productive = "T", v_call = "TRBV1*01",
    junction_aa = c("CASSF", "NA", "None"), duplicate_count = c(2, 3, 4),
    umi_count = "NA", custom_note = "NA"
  ), path, sep = "\t", quote = FALSE, row.names = FALSE)
  result <- read_airr_bulk(path, c("v_call", "junction_aa"),
    column_profile = "full", output_folder = tempfile()) |> annotations()
  expect_equal(result$junction_aa, "CASSF")
  expect_equal(result$v_call, "TRBV1*01")
  expect_equal(result$duplicate_count, 2)
  expect_true(is.numeric(result$umi_count))
  expect_true(is.na(result$umi_count))
  expect_equal(result$custom_note, "NA")
})

test_that("count normalization keeps integer-valued abundance above R integer range", {
  path <- tempfile(fileext = ".tsv")
  writeLines(c("productive\tv_call\tjunction_aa\tduplicate_count",
    "T\tTRBV1*01\tCASSF\t3000000000"), path)
  result <- read_airr_bulk(path, c("v_call", "junction_aa"),
    output_folder = tempfile()) |> annotations()
  expect_equal(result$duplicate_count, 3000000000)
  expect_equal(result$imd_n_chains, 3000000000)
})

test_that("counted bulk inputs reject missing or invalid abundance without substitution", {
  path <- tempfile(fileext = ".tsv")
  for (value in c("NA", "None", "-1", "2.5", "oops", "Inf")) {
    writeLines(c("productive\tv_call\tjunction_aa\tduplicate_count",
      paste("T", "TRBV1*01", "CASSF", value, sep = "\t")), path)
    expect_error(read_airr_bulk(path, c("v_call", "junction_aa"),
      output_folder = tempfile()), "count|Count", info = value)
  }
})

test_that("missing selected single-cell support does not silently switch count units", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    cell_id = c("cell-a", "cell-b"), locus = "IGH", productive = "T",
    v_call = "IGHV1*01", junction_aa = c("CASSF", "CASRF"),
    umi_count = c("2", "NA"), consensus_count = c(40, 90)
  ), path, row.names = FALSE)
  expect_error(read_immcantation(path, c("v_call", "junction_aa"),
    mode = "sc", output_folder = tempfile()), "Missing values.*umi_count")
})
