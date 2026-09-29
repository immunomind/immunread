test_that("10x 9.0.1 fixture maps fields and counts receptors", {
  idata <- read_fixture(read_10x, c("10x", "9.0.1", "contigs.csv.gz"))
  result <- annotations(idata)
  expect_s3_class(idata, "ImmunData")
  expect_true(all(c("cell_id", "locus", "v_call", "j_call", "junction", "junction_aa", "umi_count", "consensus_count") %in% names(result)))
  expect_true(any(!is.na(result$v_call) & result$v_call != ""))
  expect_true(any(!is.na(result$junction_aa) & result$junction_aa != ""))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 3L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 3L)
})

expected <- c("9.0.0" = 109L, "9.0.1" = 3L,
              "9.0.1-bigTRUE" = 3L, "9.0.1-camelTrue" = 346L)
for (version in names(expected)) {
  test_that(paste("10x", version, "fixture has expected receptor fields and count"), {
    idata <- read_fixture(read_10x, c("10x", version, "contigs.csv.gz"), schema = c("v_call", "j_call"))
    expect_s3_class(idata, "ImmunData")
    expect_true(all(c("cell_id", "locus", "v_call", "j_call") %in% names(annotations(idata))))
    expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), unname(expected[[version]]))
  })
}

test_that("10x reader drops unproductive contigs", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    barcode = c("cell-a", "cell-b"), chain = "TRB",
    v_gene = "TRBV1*01", cdr3 = c("CASSF", "CASRF"),
    productive = c("True", "False"), umis = c(2L, 9L)
  ), path, row.names = FALSE)

  idata <- read_10x(path, c("v_call", "junction_aa"), output_folder = tempfile())
  expect_equal(annotations(idata)$cell_id, "cell-a")
})
