# Tests for ImportRRA and BuildRRAdz functions

test_that("ImportRRA imports RRA results from directory structure", {
  temp_dir <- withr::local_tempdir()

  # Create comparison directories
  comp1_dir <- file.path(temp_dir, "d0_condA")
  comp2_dir <- file.path(temp_dir, "d0_condB")
  dir.create(comp1_dir)
  dir.create(comp2_dir)

  # Create mock RRA gene_summary output files
  rra_data1 <- data.frame(
    id = c("Trp53", "Brca1", "Myc", "Olfr123", "Gm12345", "Vmn1r50"),
    num = c(4, 5, 3, 4, 4, 3),
    pos_lfc = c(2.5, -1.8, 0.5, 1.2, -0.3, 0.8),
    pos_p_value = c(0.01, 0.03, 0.4, 0.1, 0.5, 0.2),
    pos_fdr = c(0.05, 0.1, 0.6, 0.2, 0.7, 0.3),
    neg_lfc = c(-0.5, 0.3, -0.2, -0.1, 0.1, -0.4),
    neg_p_value = c(0.5, 0.3, 0.6, 0.8, 0.9, 0.4),
    neg_fdr = c(0.7, 0.5, 0.8, 0.9, 0.95, 0.6)
  )

  rra_data2 <- data.frame(
    id = c("Trp53", "Brca1", "Myc", "Olfr456"),
    num = c(4, 5, 3, 4),
    pos_lfc = c(1.5, -2.1, 0.8, -0.5),
    pos_p_value = c(0.02, 0.01, 0.3, 0.15),
    pos_fdr = c(0.08, 0.05, 0.5, 0.25),
    neg_lfc = c(-0.3, 0.5, -0.1, 0.2),
    neg_p_value = c(0.6, 0.2, 0.7, 0.5),
    neg_fdr = c(0.8, 0.4, 0.85, 0.7)
  )

  readr::write_tsv(rra_data1, file.path(comp1_dir, "mageck_gene_summary.txt"))
  readr::write_tsv(rra_data2, file.path(comp2_dir, "mageck_gene_summary.txt"))

  # Test import
  result <- ImportRRA(temp_dir)

  # Check structure
  expect_type(result, "list")
  expect_length(result, 2)

  # Check names
  expect_true(all(c("d0_condA", "d0_condB") %in% names(result)))

  # Check filtering of Olfr, Vmn, Gm genes
  expect_false(any(grepl("^Olfr", result$d0_condA$id)))
  expect_false(any(grepl("^Vmn", result$d0_condA$id)))
  expect_false(any(grepl("^Gm\\d{4,5}", result$d0_condA$id)))

  # Check that valid genes are retained
  expect_true("Trp53" %in% result$d0_condA$id)
  expect_true("Brca1" %in% result$d0_condA$id)
})


test_that("ImportRRA handles extra_prefix parameter", {
  temp_dir <- withr::local_tempdir()
  comp_dir <- file.path(temp_dir, "d0_test")
  dir.create(comp_dir)

  rra_data <- data.frame(
    id = c("Trp53", "Brca1"),
    num = c(4, 5),
    pos_lfc = c(2.5, -1.8),
    pos_p_value = c(0.01, 0.03),
    pos_fdr = c(0.05, 0.1),
    neg_lfc = c(-0.5, 0.3),
    neg_p_value = c(0.5, 0.3),
    neg_fdr = c(0.7, 0.5)
  )

  readr::write_tsv(rra_data, file.path(comp_dir, "mageck_gene_summary.txt"))

  result <- ImportRRA(temp_dir, extra_prefix = "screen1_")

  expect_true("screen1_d0_test" %in% names(result))
})


test_that("ImportRRA converts hyphens to underscores in comparison names", {
  temp_dir <- withr::local_tempdir()
  comp_dir <- file.path(temp_dir, "day0-condA-vs-condB")
  dir.create(comp_dir)

  rra_data <- data.frame(
    id = c("Trp53"),
    num = c(4),
    pos_lfc = c(2.5),
    pos_p_value = c(0.01),
    pos_fdr = c(0.05),
    neg_lfc = c(-0.5),
    neg_p_value = c(0.5),
    neg_fdr = c(0.7)
  )

  readr::write_tsv(rra_data, file.path(comp_dir, "mageck_gene_summary.txt"))

  result <- ImportRRA(temp_dir)

  # Hyphens should be converted to underscores
  expect_true("day0_condA_vs_condB" %in% names(result))
})


test_that("BuildRRAdz builds dataframe from RRA day 0 comparisons", {
  # Create mock RRA list (simulating ImportRRA output)
  rra_list <- list(
    d0_condA = data.frame(
      id = c("Trp53", "Brca1", "Myc"),
      num = c(4, 5, 3),
      pos_lfc = c(2.5, -1.8, 0.5),
      pos_p_value = c(0.01, 0.03, 0.4),
      pos_fdr = c(0.05, 0.1, 0.6),
      neg_p_value = c(0.5, 0.2, 0.7),
      neg_fdr = c(0.7, 0.4, 0.85)
    ),
    d0_condB = data.frame(
      id = c("Trp53", "Brca1", "Myc"),
      num = c(4, 5, 3),
      pos_lfc = c(1.5, -2.1, 0.8),
      pos_p_value = c(0.02, 0.01, 0.3),
      pos_fdr = c(0.08, 0.05, 0.5),
      neg_p_value = c(0.6, 0.15, 0.5),
      neg_fdr = c(0.8, 0.3, 0.7)
    ),
    day5_condA = data.frame(
      id = c("Trp53", "Brca1"),
      num = c(4, 5),
      pos_lfc = c(3.0, -1.5),
      pos_p_value = c(0.005, 0.05),
      pos_fdr = c(0.02, 0.15),
      neg_p_value = c(0.4, 0.25),
      neg_fdr = c(0.6, 0.45)
    )
  )

  result <- BuildRRAdz(rra_list)

  # Should only include day 0 comparisons
  expect_true("id" %in% names(result))
  expect_true("num" %in% names(result))

  # Check that d0_ comparisons are included (LFC renamed with prefix)
  expect_true(any(grepl("condA_lfc", names(result))))
  expect_true(any(grepl("condB_lfc", names(result))))

  # Check minimum p-value and fdr columns
  expect_true(any(grepl("condA_pval", names(result))))
  expect_true(any(grepl("condA_fdr", names(result))))
})


test_that("BuildRRAdz requires named list", {
  unnamed_list <- list(
    data.frame(
      id = "A", num = 3, pos_lfc = 1,
      pos_p_value = 0.1, pos_fdr = 0.2,
      neg_p_value = 0.5, neg_fdr = 0.7
    ),
    data.frame(
      id = "B", num = 4, pos_lfc = 2,
      pos_p_value = 0.05, pos_fdr = 0.1,
      neg_p_value = 0.6, neg_fdr = 0.8
    )
  )

  expect_error(BuildRRAdz(unnamed_list), "rra_list must be a named list")
})


test_that("BuildRRAdz errors when no matching comparisons found", {
  rra_list <- list(
    day5_condA = data.frame(
      id = c("Trp53"),
      num = c(4),
      pos_lfc = c(3.0),
      pos_p_value = c(0.005),
      pos_fdr = c(0.02),
      neg_p_value = c(0.4),
      neg_fdr = c(0.6)
    )
  )

  expect_error(BuildRRAdz(rra_list), "No RRA objects found matching pattern")
})


test_that("BuildRRAdz respects order parameter", {
  rra_list <- list(
    d0_condA = data.frame(
      id = c("Trp53", "Brca1"),
      num = c(4, 5),
      pos_lfc = c(2.5, -1.8),
      pos_p_value = c(0.01, 0.03),
      pos_fdr = c(0.05, 0.1),
      neg_p_value = c(0.5, 0.3),
      neg_fdr = c(0.7, 0.5)
    ),
    d0_condB = data.frame(
      id = c("Trp53", "Brca1"),
      num = c(4, 5),
      pos_lfc = c(1.5, -2.1),
      pos_p_value = c(0.02, 0.01),
      pos_fdr = c(0.08, 0.05),
      neg_p_value = c(0.6, 0.2),
      neg_fdr = c(0.8, 0.4)
    )
  )

  # Test with character order
  result_ordered <- BuildRRAdz(rra_list, order = c("condB", "condA"))

  # condB columns should come before condA columns
  col_names <- names(result_ordered)
  condB_pos <- which(grepl("condB", col_names))[1]
  condA_pos <- which(grepl("condA", col_names))[1]

  expect_lt(condB_pos, condA_pos)
})


test_that("BuildRRAdz errors on invalid order parameter type", {
  rra_list <- list(
    d0_condA = data.frame(
      id = c("Trp53"),
      num = c(4),
      pos_lfc = c(2.5),
      pos_p_value = c(0.01),
      pos_fdr = c(0.05),
      neg_p_value = c(0.5),
      neg_fdr = c(0.7)
    )
  )

  expect_error(
    BuildRRAdz(rra_list, order = 123),
    "order parameter must be a character vector or factor"
  )
})


test_that("BuildRRAdz works with pattern = NULL", {
  rra_list <- list(
    exp1 = data.frame(
      id = c("Trp53", "Brca1"),
      num = c(4, 5),
      pos_lfc = c(2.5, -1.8),
      pos_p_value = c(0.01, 0.03),
      pos_fdr = c(0.05, 0.1),
      neg_p_value = c(0.5, 0.3),
      neg_fdr = c(0.7, 0.5)
    ),
    exp2 = data.frame(
      id = c("Trp53", "Brca1"),
      num = c(4, 5),
      pos_lfc = c(1.5, -2.1),
      pos_p_value = c(0.02, 0.01),
      pos_fdr = c(0.08, 0.05),
      neg_p_value = c(0.6, 0.2),
      neg_fdr = c(0.8, 0.4)
    )
  )

  # With pattern = NULL, should use all items and full names as prefixes
  result <- BuildRRAdz(rra_list, pattern = NULL)

  expect_true(any(grepl("exp1", names(result))))
  expect_true(any(grepl("exp2", names(result))))
})
