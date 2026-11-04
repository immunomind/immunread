cases <- list(
  analyzer = list(junction = TRUE, sequence = TRUE),
  analyzer_v2 = list(junction = FALSE, sequence = TRUE),
  connector = list(junction = TRUE, sequence = TRUE)
)

for (format in names(cases)) {
  test_that(paste("Adaptive", format, "fixture retains productive bulk clones"), {
    path <- fixture_path("adaptive", format, "contigs.tsv.gz")
    manifest <- data.frame(file = path, study_label = format)
    idata <- read_adaptive("<manifest>", c("v_call", "junction_aa"),
      manifest = manifest, output_folder = tempfile())
    result <- annotations(idata)

    expect_s3_class(idata, "ImmunData")
    expect_true(all(c("v_call", "d_call", "j_call", "junction_aa",
                      "duplicate_count", "productive", "study_label") %in% names(result)),
      info = format)
    expect_equal(idata |> dplyr::count() |> dplyr::pull(n), 18L,
      info = format)
    expect_equal(idata$receptors |> dplyr::count() |> dplyr::pull(n), 18L,
      info = format)
    expect_equal(idata$annotations |> dplyr::summarise(total = sum(duplicate_count)) |>
                   dplyr::pull(total), 2222, info = format)
    expect_true(all(result$productive), info = format)
    expect_true(all(result$study_label == format), info = format)
    expect_equal(result$v_call[1], "TCRBV07-09*01", info = format)
    expect_equal(result$j_call[1], "TCRBJ02-03*01", info = format)
    expect_equal(result$junction_aa[1], "CASSLARPPQWDTQYF", info = format)
    expect_identical("junction" %in% names(result), cases[[format]]$junction,
      info = format)
    expect_identical("sequence" %in% names(result), cases[[format]]$sequence,
      info = format)
    if (format == "connector") {
      expect_equal(result$junction[1],
        "TGTGCCAGCAGCTTAGCGCGTCCCCCACAGTGGGATACGCAGTATTTT")
    }
  })
}

test_that("Adaptive gene composition is not exposed as AIRR sequence origin", {
  path <- tempfile(fileext = ".tsv")
  utils::write.table(data.frame(
    v_resolved = "TRBV1*01", cdr3_amino_acid = "CASSF",
    cdr3_nucleotide = "TGTGCTTCTTCTTTT", nucleotide = "AAATGTGCTTCTTCTTTTGGG",
    templates = 12, frame_type = "In", rearrangement_type = "VDJ",
    seq_reads = 87
  ), path, sep = "\t", quote = FALSE, row.names = FALSE)
  result <- read_adaptive(path, c("v_call", "junction_aa"),
    column_profile = "full", output_folder = tempfile()) |> annotations()
  expect_false("rearrangement_type" %in% names(result))
  expect_equal(result$source_rearrangement_type, "VDJ")
  expect_equal(result$sequence, "AAATGTGCTTCTTCTTTTGGG")
  expect_equal(result$junction, "TGTGCTTCTTCTTTT")
  expect_equal(result$duplicate_count, 12)
  expect_equal(result$seq_reads, 87)
  expect_false(any(c("umi_count", "consensus_count") %in% names(result)))
})
