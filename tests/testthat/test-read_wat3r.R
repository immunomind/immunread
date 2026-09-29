test_that("wat3r 1.3 fixture maps fields and counts receptors", {
  idata <- read_fixture(read_wat3r, c("wat3r", "1.3", "contigs.csv.gz"))
  result <- annotations(idata)
  expect_s3_class(idata, "ImmunData")
  expect_true(all(c("cell_id", "locus", "v_call", "j_call", "junction", "junction_aa", "umi_count") %in% names(result)))
  expect_true(any(!is.na(result$v_call) & result$v_call != ""))
  expect_true(any(!is.na(result$junction_aa) & result$junction_aa != ""))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 1000L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 831L)
})

test_that("WAT3R requires UMI counts for each present chain", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    BC = "cell-a", TRA_CDR3 = "CAVF", TRAV = "TRAV1*01",
    TRA_CDR3_UMIcount = 2L, TRB_CDR3 = "CASSF", TRBV = "TRBV1*01"
  ), path, row.names = FALSE)
  expect_error(read_wat3r(path, c("v_call", "junction_aa"),
    output_folder = tempfile()), "TRB_CDR3_UMIcount")
})

test_that("wat3r wide records form paired TRA/TRB receptors", {
  schema <- immundata::make_receptor_schema(c("v_call", "junction_aa"), chains = c("TRA", "TRB"))
  idata <- read_fixture(read_wat3r, c("wat3r", "1.3", "contigs.csv.gz"), schema = schema)
  expect_setequal(unique(annotations(idata)$locus), c("TRA", "TRB"))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 808L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 388L)
})
