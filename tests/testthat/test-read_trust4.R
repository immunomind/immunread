test_that("trust4 1.1.7 fixture maps fields and counts receptors", {
  idata <- read_fixture(read_trust4, c("trust4", "1.1.7", "contigs.csv.gz"))
  result <- annotations(idata)
  expect_s3_class(idata, "ImmunData")
  expect_true(all(c("cell_id", "locus", "v_call", "j_call", "junction", "junction_aa") %in% names(result)))
  expect_true(any(!is.na(result$v_call) & result$v_call != ""))
  expect_true(any(!is.na(result$junction_aa) & result$junction_aa != ""))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 1000L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 845L)
})

test_that("TRUST4 uses its AIRR export count convention and keeps the estimate", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    X.barcode = "cell-a",
    chain1 = "TRBV1*01,*,TRBJ1*01,TRBC1,TGTGCTTCTTCTTTT,CASSF,7.90,contig-a",
    chain2 = "*", secondary_chain1 = "*", secondary_chain2 = "*"
  ), path, row.names = FALSE)
  result <- read_trust4(path, c("v_call", "junction_aa"),
    column_profile = "full", output_folder = tempfile()) |> annotations()
  expect_equal(result$consensus_count, 7)
  expect_equal(result$source_read_count, 7.9)
  expect_equal(result$junction, "TGTGCTTCTTCTTTT")
  expect_equal(result$junction_aa, "CASSF")
  expect_true(is.na(result$d_call))
})
