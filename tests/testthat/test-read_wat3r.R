test_that("wat3r 1.3 fixture maps fields and counts receptors", {
  idata <- read_fixture(read_wat3r, c("wat3r", "1.3", "contigs.csv.gz"))
  result <- annotations(idata)
  expect_s3_class(idata, "ImmunData")
  expect_true(all(c("cell_id", "locus", "v_call", "j_call", "cdr3", "junction_aa", "umi_count") %in% names(result)))
  expect_true(any(!is.na(result$v_call) & result$v_call != ""))
  expect_true(any(!is.na(result$junction_aa) & result$junction_aa != ""))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 1000L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 831L)
  expect_false("junction" %in% names(result))
})

test_that("WAT3R maps CDR3 sequence and read and UMI support separately", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    BC = c("cell-tra", "cell-trb", "cell-tra2"),
    TRAV = c("TRAV1*01", NA, NA), TRAJ = c("TRAJ1*01", NA, NA),
    TRA_CDR3 = c("CAVF", NA, NA), TRA_CDR3nuc = c("GCTGTA", NA, NA),
    TRA_nReads = c(11L, NA, NA), TRA_CDR3_UMIcount = c(3L, NA, NA),
    TRBV = c(NA, "TRBV1*01", NA), TRBJ = c(NA, "TRBJ1*01", NA),
    TRB_CDR3 = c(NA, "CASSF", NA), TRB_CDR3nuc = c(NA, "GCCAGCAGT", NA),
    TRB_nReads = c(NA, 22L, NA), TRB_CDR3_UMIcount = c(NA, 5L, NA),
    `TRAV.2` = c(NA, NA, "TRAV2*01"), `TRAJ.2` = c(NA, NA, "TRAJ2*01"),
    TRA.2_CDR3 = c(NA, NA, "CAVVF"), TRA.2_CDR3nuc = c(NA, NA, "GCTGTAGTT"),
    TRA.2_nReads = c(NA, NA, 33L), TRA.2_CDR3_UMIcount = c(NA, NA, 7L)
  ), path, row.names = FALSE)
  idata <- read_wat3r(path, c("v_call", "junction_aa", "cdr3"),
    output_folder = tempfile())
  result <- annotations(idata)
  expect_equal(result |> dplyr::arrange(cell_id) |> dplyr::pull(cdr3),
    c("GCTGTA", "GCTGTAGTT", "GCCAGCAGT"))
  expect_equal(result |> dplyr::arrange(cell_id) |> dplyr::pull(junction_aa),
    c("CAVF", "CAVVF", "CASSF"))
  expect_false("junction" %in% names(result))
  expect_equal(result |> dplyr::arrange(cell_id) |> dplyr::pull(consensus_count),
    c(11, 33, 22))
  expect_equal(result |> dplyr::arrange(cell_id) |> dplyr::pull(umi_count),
    c(3, 7, 5))

  no_reads_path <- tempfile(fileext = ".csv")
  no_reads <- data.frame(
    BC = "cell-no-reads", TRAV = "TRAV1*01", TRAJ = "TRAJ1*01",
    TRA_CDR3 = "CAVF", TRA_CDR3nuc = "GCTGTA",
    TRA_CDR3_UMIcount = 3L
  )
  utils::write.csv(no_reads, no_reads_path, row.names = FALSE)
  no_reads_idata <- read_wat3r(no_reads_path, c("v_call", "junction_aa", "cdr3"),
    output_folder = tempfile())
  expect_false("consensus_count" %in% names(annotations(no_reads_idata)))
})

test_that("wat3r wide records form paired TRA/TRB receptors", {
  schema <- immundata::make_receptor_schema(c("v_call", "junction_aa"), chains = c("TRA", "TRB"))
  idata <- read_fixture(read_wat3r, c("wat3r", "1.3", "contigs.csv.gz"), schema = schema)
  expect_setequal(unique(annotations(idata)$locus), c("TRA", "TRB"))
  expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 808L)
  expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 388L)
})
