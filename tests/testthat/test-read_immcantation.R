test_that("Immcantation single-cell fixture preserves UMI and read support", {
  path <- fixture_path("immcantation", "single-cell.tsv.gz")
  input <- utils::read.delim(path, check.names = FALSE)
  expect_equal(nrow(input), 80L)
  expect_equal(length(unique(input$cell_id)), 40L)
  expect_false("duplicate_count" %in% names(input))

  schema <- immundata::make_receptor_schema(
    c("v_call", "junction_aa"), chains = "IGH"
  )
  idata <- read_immcantation_sc(path, schema, output_folder = tempfile())
  result <- annotations(idata)
  expect_s3_class(idata, "ImmunData")
  expect_equal(nrow(result), 38L)
  expect_equal(length(unique(result$cell_id)), 38L)
  expect_false(anyNA(result$junction_aa))
  expect_false(any(result$junction_aa %in% c("", "NA", "None")))

  # The input contains two heavy chains in this cell, supported by 46 and 31 UMIs.
  heavy <- result[result$cell_id == "TCGTACCGTGACGCCT-1" & result$locus == "IGH", ]
  expect_equal(heavy$junction_aa, "CVKDRSSRQLKAVGTYFDCW")
  expect_equal(heavy$umi_count, 46)
  expect_equal(heavy$consensus_count, 2529)
})

test_that("Immcantation bulk fixture preserves sequence copy abundance", {
  idata <- read_fixture(
    read_immcantation_bulk, c("immcantation", "bulk.tsv.gz"),
    schema = c("v_call", "junction"),
    repertoire_schema = "sample_id"
  )
  result <- annotations(idata)
  expect_s3_class(idata, "ImmunData")
  expect_equal(nrow(result), 1780L)
  expect_setequal(unique(result$sample_id), c("-1h", "+7d"))
  expect_equal(sum(result$duplicate_count), 4148)
  expect_equal(result$imd_n_chains, result$duplicate_count)
  expect_false(any(c("umi_count", "consensus_count", "junction_aa") %in% names(result)))

  # A known high-copy observation protects against silently treating rows as counts.
  high_copy <- result[result$sequence_id == "GN5SHBT01A55XE", ]
  expect_equal(high_copy$duplicate_count, 464)
  expect_equal(high_copy$junction,
    "TGTAGTAGAGATCTCGCGGTTATATCCACAATAGCTGGTACTAACTGGTTCGACCCCAGG")
})

test_that("Immcantation single-cell requires a real UMI or read count", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    cell_id = "cell-a", locus = "TRB", productive = "T",
    v_call = "TRBV1*01", junction_aa = "CASSF"
  ), path, row.names = FALSE)
  expect_error(read_immcantation_sc(path, c("v_call", "junction_aa"),
    output_folder = tempfile()), "requires.*umi_col")
})
