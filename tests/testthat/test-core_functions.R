test_that("parse_barcode_policy works correctly", {
  # Valid policy
  policy <- parse_barcode_policy("PREFIX:CACCG@6")
  expect_equal(policy$flanking_seq, "CACCG")
  expect_equal(policy$position_offset, 6)

  # Invalid policy - wrong format
  expect_error(
    parse_barcode_policy("INVALID:FORMAT"),
    "Invalid barcode policy format"
  )

  # Invalid policy - non-DNA sequence
  expect_error(
    parse_barcode_policy("PREFIX:ACCTGN@6"),
    "Invalid barcode policy format"
  )
})


test_that("extract_barcodes_from_read extracts correct barcodes", {
  # Mock barcode policy
  policy <- list(flanking_seq = "CACCG", position_offset = 6)

  # Perfect read with flanking sequence and barcode
  # Sample barcode: TTGAACCG, Flanking: CACCG at position 9, then 20bp guide
  sequence <- "TTGAACCGCACCGAAAAAAAAAAAATGCATTCTNNNNNNNN"
  quality <- paste(rep("I", nchar(sequence)), collapse = "")

  result <- extract_barcodes_from_read(
    sequence,
    quality,
    "read_001",
    policy,
    construct_barcode_length = 20
  )

  expect_equal(result$read_id, "read_001")
  expect_equal(result$sample_barcode, "TTGAACCG")
  expect_equal(result$construct_barcode, "AAAAAAAAAAAATGCATTCT")
  expect_equal(result$extraction_status, "success")
})


test_that("create_count_matrix handles construct level", {
  # Create minimal test data with required columns including condition
  test_data <- data.table::data.table(
    barcode = c("AAAAAAAAAAAATGCATTCT", "AAAAAAAAAAAATGCATTCT", "AAAAAAAAATAAGCTCACCC"),
    construct_id = c("BRDN0001468998", "BRDN0001468998", "BRDN0001474573"),
    gene_symbol = c("Gene1", "Gene1", "Gene2"),
    gene_id = c("1001", "1001", "1002"),
    sample_barcode = c("TTGAACCG", "TTGAACCG", "AATCCAGC"),
    condition = c("Sample1", "Sample1", "Sample2")
  )

  result <- create_count_matrix(test_data, level = "construct", quiet = TRUE)

  expect_type(result, "list")
  expect_true(is.matrix(result$matrix))
  expect_equal(nrow(result$matrix), 2)  # 2 unique constructs
  expect_equal(ncol(result$matrix), 2)  # 2 unique samples

  # Check that counts are correct
  expect_equal(as.integer(result$matrix["BRDN0001468998", "Sample1"]), 2L)
  expect_equal(as.integer(result$matrix["BRDN0001474573", "Sample2"]), 1L)
})


test_that("lognormalize_counts applies poolq normalization", {
  # Create a simple count matrix
  count_matrix <- matrix(
    c(146, 0, 228, 215, 5, 231),
    nrow = 2,
    ncol = 3,
    byrow = FALSE,
    dimnames = list(
      c("construct1", "construct2"),
      c("sample1", "sample2", "sample3")
    )
  )

  # Apply normalization
  lognorm <- lognormalize_counts(count_matrix, method = "poolq", quiet = TRUE)

  # Check properties
  expect_type(lognorm, "double")
  expect_equal(dim(lognorm), dim(count_matrix))

  # Check that values are reasonable (should be mostly positive/negative depending on depth)
  expect_true(all(is.numeric(lognorm)))
  expect_true(all(is.finite(lognorm)))
})


test_that("calculate_qc_metrics handles valid input", {
  # Create minimal test data
  extracted <- data.table::data.table(
    read_id = c("r1", "r2", "r3", "r4"),
    construct_barcode = c("barcode1", "barcode2", "barcode3", NA),
    sample_barcode = c("TTGAACCG", "TTGAACCG", "AATCCAGC", NA),
    extraction_status = c("success", "success", "success", "no_construct_flank")
  )

  matched <- extracted[!is.na(construct_barcode)]
  matched[, construct_id := c("BRDN0001468998", "BRDN0001474573", "BRDN0001438906")]

  aggregated <- data.table::copy(matched)
  aggregated[, `:=`(
    condition = c("Day0", "Day0", "UT_r1"),
    gene_id = c("1001", "1002", "1003"),
    gene_symbol = c("Gene1", "Gene2", "Gene3")
  )]

  qc <- calculate_qc_metrics(extracted, matched, aggregated, quiet = TRUE)

  expect_type(qc, "list")
  expect_type(qc$overall_metrics, "list")
  expect_type(qc$per_sample_metrics, "list")

  # Check overall metrics
  expect_true(nrow(qc$overall_metrics) > 0)
  expect_true("metric" %in% names(qc$overall_metrics))

  # Check per-sample metrics
  expect_true(nrow(qc$per_sample_metrics) > 0)
})
