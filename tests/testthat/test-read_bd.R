test_that("bd 2.3 fixture maps fields and counts receptors", {
  idata <- read_fixture(read_bd, c("bd", "2.3", "contigs.csv.gz"))
  result <- annotations(idata)
  expect_s3_class(idata, "ImmunData")
  expect_true(all(c("cell_id", "locus", "v_call", "j_call", "junction", "junction_aa", "umi_count", "consensus_count") %in% names(result)))
  expect_true(any(!is.na(result$v_call) & result$v_call != ""))
  expect_true(any(!is.na(result$junction_aa) & result$junction_aa != ""))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 108L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 106L)
})

test_that("BD normalizes old and modern UMI count fields", {
  write_input <- function(data) {
    path <- tempfile(fileext = ".csv")
    utils::write.csv(data, path, row.names = FALSE)
    path
  }
  base <- data.frame(
    cell_id = c("cell-a", "cell-b"), locus = "TRB",
    v_call = "TRBV1*01", junction_aa = c("CASSF", "CASSG"),
    productive = "T", consensus_count = c(20L, 30L)
  )

  old <- base
  old$duplicate_count <- c(2L, 8L)
  old_result <- annotations(read_bd(write_input(old), "junction_aa",
    column_profile = "full", output_folder = tempfile())) |> dplyr::arrange(cell_id)
  expect_true(all(c("umi_count", "consensus_count", "source_duplicate_count") %in% names(old_result)))
  expect_equal(old_result$umi_count, c(2L, 8L))
  expect_equal(old_result$consensus_count, c(20L, 30L))

  modern <- base
  modern$umi_count <- c(3L, 7L)
  modern$duplicate_count <- c(11L, 13L)
  modern_result <- annotations(read_bd(write_input(modern), "junction_aa",
    column_profile = "full", output_folder = tempfile())) |> dplyr::arrange(cell_id)
  expect_equal(modern_result$umi_count, c(3L, 7L))
  expect_equal(modern_result$duplicate_count, c(11L, 13L))
})
