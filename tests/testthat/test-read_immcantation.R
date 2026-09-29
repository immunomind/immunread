test_that("immcantation 4.6.0 fixture maps fields and counts receptors", {
  idata <- read_fixture(read_immcantation, c("immcantation", "4.6.0", "contigs.csv.gz"), mode = "sc")
  result <- annotations(idata)
  expect_s3_class(idata, "ImmunData")
  expect_true(all(c("cell_id", "locus", "v_call", "j_call", "junction", "junction_aa", "duplicate_count") %in% names(result)))
  expect_true(any(!is.na(result$v_call) & result$v_call != ""))
  expect_true(any(!is.na(result$junction_aa) & result$junction_aa != ""))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 174L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 70L)
})

test_that("Immcantation single-cell requires a real UMI or read count", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    cell_id = "cell-a", locus = "TRB", productive = "T",
    v_call = "TRBV1*01", junction_aa = "CASSF"
  ), path, row.names = FALSE)
  expect_error(read_immcantation(path, c("v_call", "junction_aa"),
    mode = "sc", output_folder = tempfile()), "requires.*umi_col")
})
