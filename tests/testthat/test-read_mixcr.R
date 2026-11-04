test_that("mixcr 4.7.0 fixture maps fields and counts receptors", {
  idata <- read_fixture(read_mixcr, c("mixcr", "4.7.0", "contigs.csv.gz"))
  result <- annotations(idata)
  expect_s3_class(idata, "ImmunData")
  expect_true(all(c("cell_id", "locus", "v_call", "j_call", "junction", "junction_aa", "umi_count", "consensus_count") %in% names(result)))
  expect_false("duplicate_count" %in% names(result))
  expect_true(any(!is.na(result$v_call) & result$v_call != ""))
  expect_true(any(!is.na(result$junction_aa) & result$junction_aa != ""))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 194L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 194L)
})

test_that("mixcr readCount maps to consensus_count and stays in full profile", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    cellGroup = "cell-a", topChains = "TRA", allVHitsWithScore = "TRAV1*01",
    aaSeqCDR3 = "CAVF",
    uniqueMoleculeCount = 3L, readCount = 9L
  ), path, row.names = FALSE)
  idata <- read_mixcr(path, c("v_call", "junction_aa"),
    column_profile = "full", output_folder = tempfile())
  result <- annotations(idata)
  expect_equal(result$umi_count, 3L)
  expect_equal(result$consensus_count, 9L)
  expect_equal(result$readCount, 9L)
  expect_false("duplicate_count" %in% names(result))
})
