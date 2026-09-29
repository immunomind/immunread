test_that("dandelion 0.5.5 fixture maps fields and counts receptors", {
  idata <- read_fixture(read_dandelion, c("dandelion", "0.5.5", "contigs.csv.gz"))
  result <- annotations(idata)
  expect_s3_class(idata, "ImmunData")
  expect_true(all(c("cell_id", "locus", "v_call", "j_call", "junction", "junction_aa", "duplicate_count") %in% names(result)))
  expect_true(any(!is.na(result$v_call) & result$v_call != ""))
  expect_true(any(!is.na(result$junction_aa) & result$junction_aa != ""))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 121L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 121L)
})

test_that("Dandelion single-cell requires a real UMI or read count", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    cell_id = "cell-a", locus = "TRB", productive = "T",
    v_call = "TRBV1*01", junction_aa = "CASSF"
  ), path, row.names = FALSE)
  expect_error(read_dandelion(path, c("v_call", "junction_aa"),
    output_folder = tempfile()), "requires.*umi_col")
})
