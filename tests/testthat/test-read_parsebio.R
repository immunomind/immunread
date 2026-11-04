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
  ), path, row.names = FALSE, na = "")
  expect_error(read_parsebio(path, c("v_call", "junction_aa"),
    output_folder = tempfile()), "TRB_transcript_count")
})

test_that("ParseBio maps transcript and read counts for every chain", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    Barcode = paste0("cell-", 1:4),
    TRA_cdr3_aa = c("CAVF", NA, NA, NA), TRA_V = c("TRAV1*01", NA, NA, NA),
    TRA_transcript_count = c(2, NA, NA, NA), TRA_read_count = c(7, NA, NA, NA),
    TRB_cdr3_aa = c(NA, "CASSF", NA, NA), TRB_V = c(NA, "TRBV1*01", NA, NA),
    TRB_transcript_count = c(NA, 3, NA, NA), TRB_read_count = c(NA, 11, NA, NA),
    secondary_TRA_cdr3_aa = c(NA, NA, "CAVX", NA),
    secondary_TRA_V = c(NA, NA, "TRAV2*01", NA),
    secondary_TRA_transcript_count = c(NA, NA, 4, NA),
    secondary_TRA_read_count = c(NA, NA, 13, NA),
    secondary_TRB_cdr3_aa = c(NA, NA, NA, "CASXX"),
    secondary_TRB_V = c(NA, NA, NA, "TRBV2*01"),
    secondary_TRB_transcript_count = c(NA, NA, NA, 5),
    secondary_TRB_read_count = c(NA, NA, NA, 17)
  ), path, row.names = FALSE, na = "")
  idata <- read_parsebio(path, c("v_call", "junction_aa"),
    output_folder = tempfile())
  result <- annotations(idata)
  expect_setequal(result$cell_id, paste0("cell-", 1:4))
  expect_setequal(result$junction_aa, c("CAVF", "CASSF", "CAVX", "CASXX"))
  expect_setequal(result$umi_count, c(2, 3, 4, 5))
  expect_setequal(result$consensus_count, c(7, 11, 13, 17))
})

test_that("ParseBio allows missing optional read counts", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    Barcode = c("cell-a", "cell-b"),
    TRA_cdr3_aa = c("CAVF", NA), TRA_V = c("TRAV1*01", NA),
    TRA_transcript_count = c(2, NA), TRA_read_count = c(7, NA),
    TRB_cdr3_aa = c(NA, "CASSF"), TRB_V = c(NA, "TRBV1*01"),
    TRB_transcript_count = c(NA, 3), TRB_read_count = c(NA, NA)
  ), path, row.names = FALSE)
  idata <- read_parsebio(path, c("v_call", "junction_aa"),
    output_folder = tempfile())
  result <- annotations(idata)
  expect_equal(result$consensus_count[result$cell_id == "cell-a"], 7)
  expect_type(result$consensus_count, "double")
  expect_true(is.na(result$consensus_count[result$cell_id == "cell-b"]))
})

test_that("parsebio wide records form paired TRA/TRB receptors", {
  schema <- immundata::make_receptor_schema(c("v_call", "junction_aa"), chains = c("TRA", "TRB"))
  idata <- read_fixture(read_parsebio, c("parsebio", "1.6.1", "contigs.csv.gz"), schema = schema)
  expect_setequal(unique(annotations(idata)$locus), c("TRA", "TRB"))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 1144L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 565L)
})
