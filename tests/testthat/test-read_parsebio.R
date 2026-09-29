test_that("parsebio 1.6.1 fixture maps fields and counts receptors", {
  idata <- read_fixture(read_parsebio, c("parsebio", "1.6.1", "contigs.csv.gz"))
  result <- annotations(idata)
  expect_s3_class(idata, "ImmunData")
  expect_true(all(c("cell_id", "locus", "v_call", "j_call", "junction_aa", "umi_count") %in% names(result)))
  expect_true(any(!is.na(result$v_call) & result$v_call != ""))
  expect_true(any(!is.na(result$junction_aa) & result$junction_aa != ""))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 1000L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 975L)
})

test_that("ParseBio requires transcript counts for each present chain", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    Barcode = "cell-a", TRA_cdr3_aa = "CAVF", TRA_V = "TRAV1*01",
    TRA_transcript_count = 2L, TRB_cdr3_aa = "CASSF", TRB_V = "TRBV1*01"
  ), path, row.names = FALSE)
  expect_error(read_parsebio(path, c("v_call", "junction_aa"),
    output_folder = tempfile()), "TRB_transcript_count")
})

test_that("parsebio wide records form paired TRA/TRB receptors", {
  schema <- immundata::make_receptor_schema(c("v_call", "junction_aa"), chains = c("TRA", "TRB"))
  idata <- read_fixture(read_parsebio, c("parsebio", "1.6.1", "contigs.csv.gz"), schema = schema)
  expect_setequal(unique(annotations(idata)$locus), c("TRA", "TRB"))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 1144L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 565L)
})
