# Tests for ImportMLE function

test_that("ImportMLE imports MLE results from directory structure", {
  temp_dir <- withr::local_tempdir()

  # Create MLE result directories with default prefix

comp1_dir <- file.path(temp_dir, "mageck_MLE_condA")
  comp2_dir <- file.path(temp_dir, "mageck_MLE_condB")
  dir.create(comp1_dir)
  dir.create(comp2_dir)

  # Create mock MLE gene_summary output files
  mle_data1 <- data.frame(
    Gene = c("Trp53", "Brca1", "Myc", "Olfr123", "Gm12345", "Vmn1r50"),
    sgRNA = c(4, 5, 3, 4, 4, 3),
    beta = c(2.5, -1.8, 0.5, 1.2, -0.3, 0.8),
    z = c(3.1, -2.5, 0.8, 1.5, -0.4, 1.0),
    `p-value` = c(0.01, 0.03, 0.4, 0.1, 0.5, 0.2),
    fdr = c(0.05, 0.1, 0.6, 0.2, 0.7, 0.3),
    check.names = FALSE
  )

  mle_data2 <- data.frame(
    Gene = c("Trp53", "Brca1", "Myc", "Olfr456"),
    sgRNA = c(4, 5, 3, 4),
    beta = c(1.5, -2.1, 0.8, -0.5),
    z = c(2.0, -3.0, 1.2, -0.7),
    `p-value` = c(0.02, 0.01, 0.3, 0.15),
    fdr = c(0.08, 0.05, 0.5, 0.25),
    check.names = FALSE
  )

  readr::write_tsv(mle_data1, file.path(comp1_dir, "mageck_gene_summary.txt"))
  readr::write_tsv(mle_data2, file.path(comp2_dir, "mageck_gene_summary.txt"))

  # Test import (suppress messages about default prefix/dir)
  result <- suppressMessages(ImportMLE(temp_dir))

  # Check structure
  expect_type(result, "list")
  expect_length(result, 2)

  # Check names (prefix stripped)
  expect_true(all(c("condA", "condB") %in% names(result)))

  # Check filtering of Olfr, Vmn, Gm genes
  expect_false(any(grepl("^Olfr", result$condA$Gene)))
  expect_false(any(grepl("^Vmn", result$condA$Gene)))
  expect_false(any(grepl("^Gm\\d{4,5}", result$condA$Gene)))

  # Check that valid genes are retained
  expect_true("Trp53" %in% result$condA$Gene)
  expect_true("Brca1" %in% result$condA$Gene)
})


test_that("ImportMLE handles extra_prefix parameter", {
  temp_dir <- withr::local_tempdir()
  comp_dir <- file.path(temp_dir, "mageck_MLE_test")
  dir.create(comp_dir)

  mle_data <- data.frame(
    Gene = c("Trp53", "Brca1"),
    sgRNA = c(4, 5),
    beta = c(2.5, -1.8),
    z = c(3.1, -2.5),
    `p-value` = c(0.01, 0.03),
    fdr = c(0.05, 0.1),
    check.names = FALSE
  )

  readr::write_tsv(mle_data, file.path(comp_dir, "mageck_gene_summary.txt"))

  result <- suppressMessages(ImportMLE(temp_dir, extra_prefix = "screen1_"))

  expect_true("screen1_test" %in% names(result))
})


test_that("ImportMLE handles custom prefix parameter", {
  temp_dir <- withr::local_tempdir()
  comp_dir <- file.path(temp_dir, "custom_prefix_condA")
  dir.create(comp_dir)

  mle_data <- data.frame(
    Gene = c("Trp53", "Brca1"),
    sgRNA = c(4, 5),
    beta = c(2.5, -1.8),
    z = c(3.1, -2.5),
    `p-value` = c(0.01, 0.03),
    fdr = c(0.05, 0.1),
    check.names = FALSE
  )

  readr::write_tsv(mle_data, file.path(comp_dir, "mageck_gene_summary.txt"))

  # Use custom prefix
  result <- suppressMessages(ImportMLE(temp_dir, prefix = "custom_prefix_"))

  expect_true("condA" %in% names(result))
})


test_that("ImportMLE shows message when using default directory", {
  temp_dir <- withr::local_tempdir()
  comp_dir <- file.path(temp_dir, "mageck_MLE_test")
  dir.create(comp_dir)

  mle_data <- data.frame(
    Gene = c("Trp53"),
    sgRNA = c(4),
    beta = c(2.5),
    z = c(3.1),
    `p-value` = c(0.01),
    fdr = c(0.05),
    check.names = FALSE
  )

  readr::write_tsv(mle_data, file.path(comp_dir, "mageck_gene_summary.txt"))

  # Using current directory should show message
  withr::with_dir(temp_dir, {
    expect_message(ImportMLE(), "Searching in the top layer of working directory")
  })
})


test_that("ImportMLE shows message when using default prefix", {
  temp_dir <- withr::local_tempdir()
  comp_dir <- file.path(temp_dir, "mageck_MLE_test")
  dir.create(comp_dir)

  mle_data <- data.frame(
    Gene = c("Trp53"),
    sgRNA = c(4),
    beta = c(2.5),
    z = c(3.1),
    `p-value` = c(0.01),
    fdr = c(0.05),
    check.names = FALSE
  )

  readr::write_tsv(mle_data, file.path(comp_dir, "mageck_gene_summary.txt"))

  expect_message(ImportMLE(temp_dir), "MLE results folders will be identified with the prefix")
})


test_that("ImportMLE cleans pipe and hyphen characters from column names", {
  temp_dir <- withr::local_tempdir()
  comp_dir <- file.path(temp_dir, "mageck_MLE_test")
  dir.create(comp_dir)

  # Create data with problematic column names
  mle_data <- data.frame(
    Gene = c("Trp53", "Brca1"),
    sgRNA = c(4, 5),
    `beta|treatment` = c(2.5, -1.8),
    `z-score` = c(3.1, -2.5),
    `p-value|adj` = c(0.01, 0.03),
    fdr = c(0.05, 0.1),
    check.names = FALSE
  )

  readr::write_tsv(mle_data, file.path(comp_dir, "mageck_gene_summary.txt"))

  result <- suppressMessages(ImportMLE(temp_dir))

  # Check that pipe and hyphen characters are replaced with underscores
  col_names <- names(result$test)
  expect_false(any(grepl("\\|", col_names)))
  expect_false(any(grepl("-", col_names)))
})


test_that("ImportMLE returns empty list when no matching directories found", {
  temp_dir <- withr::local_tempdir()

  # Create directory that doesn't match the prefix
  other_dir <- file.path(temp_dir, "other_results")
  dir.create(other_dir)

  # Should return empty list
  result <- suppressMessages(ImportMLE(temp_dir))

  expect_type(result, "list")
  expect_length(result, 0)
})
