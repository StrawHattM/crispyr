# Tests for ImportDrugZ and BuildDrugZdz functions

test_that("ImportDrugZ imports DrugZ results from directory structure", {
  # Create temporary directory structure mimicking DrugZ output

  temp_dir <- withr::local_tempdir()

  # Create comparison directories
  comp1_dir <- file.path(temp_dir, "d0_condA")
  comp2_dir <- file.path(temp_dir, "d0_condB")
  dir.create(comp1_dir)
  dir.create(comp2_dir)

  # Create mock DrugZ output files
  drugz_data1 <- data.frame(
    GENE = c("Trp53", "Brca1", "Myc", "Olfr123", "Gm12345", "Vmn1r50"),
    numObs = c(4, 5, 3, 4, 4, 3),
    normZ = c(2.5, -1.8, 0.5, 1.2, -0.3, 0.8),
    pval_synth = c(0.01, 0.03, 0.4, 0.1, 0.5, 0.2),
    fdr_synth = c(0.05, 0.1, 0.6, 0.2, 0.7, 0.3)
  )

  drugz_data2 <- data.frame(
    GENE = c("Trp53", "Brca1", "Myc", "Olfr456"),
    numObs = c(4, 5, 3, 4),
    normZ = c(1.5, -2.1, 0.8, -0.5),
    pval_synth = c(0.02, 0.01, 0.3, 0.15),
    fdr_synth = c(0.08, 0.05, 0.5, 0.25)
  )

  readr::write_tsv(drugz_data1, file.path(comp1_dir, "drugz_output.txt"))
  readr::write_tsv(drugz_data2, file.path(comp2_dir, "drugz_output.txt"))

  # Test import
  result <- ImportDrugZ(temp_dir)

  # Check structure

  expect_type(result, "list")
  expect_length(result, 2)

  # Check names (hyphens converted to underscores)
  expect_true(all(c("d0_condA", "d0_condB") %in% names(result)))

  # Check filtering of Olfr, Vmn, Gm genes
  expect_false(any(grepl("^Olfr", result$d0_condA$GENE)))
  expect_false(any(grepl("^Vmn", result$d0_condA$GENE)))
  expect_false(any(grepl("^Gm\\d{4,5}", result$d0_condA$GENE)))

  # Check that valid genes are retained
  expect_true("Trp53" %in% result$d0_condA$GENE)
  expect_true("Brca1" %in% result$d0_condA$GENE)
})


test_that("ImportDrugZ handles extra_prefix parameter", {
  temp_dir <- withr::local_tempdir()
  comp_dir <- file.path(temp_dir, "d0_test")

  dir.create(comp_dir)

  drugz_data <- data.frame(
    GENE = c("Trp53", "Brca1"),
    numObs = c(4, 5),
    normZ = c(2.5, -1.8),
    pval_synth = c(0.01, 0.03),
    fdr_synth = c(0.05, 0.1)
  )

  readr::write_tsv(drugz_data, file.path(comp_dir, "drugz_output.txt"))

  result <- ImportDrugZ(temp_dir, extra_prefix = "exp1_")

  expect_true("exp1_d0_test" %in% names(result))
})


test_that("ImportDrugZ converts hyphens to underscores in comparison names", {
  temp_dir <- withr::local_tempdir()
  comp_dir <- file.path(temp_dir, "day0-condA-vs-condB")
  dir.create(comp_dir)

  drugz_data <- data.frame(
    GENE = c("Trp53"),
    numObs = c(4),
    normZ = c(2.5),
    pval_synth = c(0.01),
    fdr_synth = c(0.05)
  )

  readr::write_tsv(drugz_data, file.path(comp_dir, "drugz_output.txt"))

  result <- ImportDrugZ(temp_dir)

  # Hyphens should be converted to underscores
  expect_true("day0_condA_vs_condB" %in% names(result))
})


test_that("BuildDrugZdz builds dataframe from DrugZ day 0 comparisons", {
  # Create mock DrugZ list (simulating ImportDrugZ output)
  drugz_list <- list(
    d0_condA = data.frame(
      GENE = c("Trp53", "Brca1", "Myc"),
      numObs = c(4, 5, 3),
      normZ = c(2.5, -1.8, 0.5),
      pval_synth = c(0.01, 0.03, 0.4),
      fdr_synth = c(0.05, 0.1, 0.6)
    ),
    d0_condB = data.frame(
      GENE = c("Trp53", "Brca1", "Myc"),
      numObs = c(4, 5, 3),
      normZ = c(1.5, -2.1, 0.8),
      pval_synth = c(0.02, 0.01, 0.3),
      fdr_synth = c(0.08, 0.05, 0.5)
    ),
    day5_condA = data.frame(
      GENE = c("Trp53", "Brca1"),
      numObs = c(4, 5),
      normZ = c(3.0, -1.5),
      pval_synth = c(0.005, 0.05),
      fdr_synth = c(0.02, 0.15)
    )
  )

  result <- BuildDrugZdz(drugz_list)

  # Should only include day 0 comparisons
  expect_true("GENE" %in% names(result))
  expect_true("numObs" %in% names(result))

  # Check that d0_ comparisons are included
  expect_true(any(grepl("condA_normZ", names(result))))
  expect_true(any(grepl("condB_normZ", names(result))))

  # Check minimum p-value calculation
  expect_true(any(grepl("condA_pval", names(result))))
  expect_true(any(grepl("condA_fdr", names(result))))
})


test_that("BuildDrugZdz requires named list", {
  unnamed_list <- list(
    data.frame(GENE = "A", numObs = 3, normZ = 1, pval_synth = 0.1, fdr_synth = 0.2),
    data.frame(GENE = "B", numObs = 4, normZ = 2, pval_synth = 0.05, fdr_synth = 0.1)
  )

  expect_error(BuildDrugZdz(unnamed_list), "drugz_list must be a named list")
})


test_that("BuildDrugZdz errors when no matching comparisons found", {
  drugz_list <- list(
    day5_condA = data.frame(
      GENE = c("Trp53"),
      numObs = c(4),
      normZ = c(3.0),
      pval_synth = c(0.005),
      fdr_synth = c(0.02)
    )
  )

  expect_error(BuildDrugZdz(drugz_list), "No DrugZ objects found matching pattern")
})


test_that("BuildDrugZdz respects order parameter", {
  drugz_list <- list(
    d0_condA = data.frame(
      GENE = c("Trp53", "Brca1"),
      numObs = c(4, 5),
      normZ = c(2.5, -1.8),
      pval_synth = c(0.01, 0.03),
      fdr_synth = c(0.05, 0.1)
    ),
    d0_condB = data.frame(
      GENE = c("Trp53", "Brca1"),
      numObs = c(4, 5),
      normZ = c(1.5, -2.1),
      pval_synth = c(0.02, 0.01),
      fdr_synth = c(0.08, 0.05)
    )
  )

  # Test with character order
  result_ordered <- BuildDrugZdz(drugz_list, order = c("condB", "condA"))

  # condB columns should come before condA columns
  col_names <- names(result_ordered)
  condB_pos <- which(grepl("condB", col_names))[1]
  condA_pos <- which(grepl("condA", col_names))[1]

  expect_lt(condB_pos, condA_pos)
})


test_that("BuildDrugZdz errors on invalid order parameter type", {
  drugz_list <- list(
    d0_condA = data.frame(
      GENE = c("Trp53"),
      numObs = c(4),
      normZ = c(2.5),
      pval_synth = c(0.01),
      fdr_synth = c(0.05)
    )
  )

  expect_error(
    BuildDrugZdz(drugz_list, order = 123),
    "order parameter must be a character vector or factor"
  )
})


test_that("BuildDrugZdz works with pattern = NULL",
          {
  drugz_list <- list(
    exp1 = data.frame(
      GENE = c("Trp53", "Brca1"),
      numObs = c(4, 5),
      normZ = c(2.5, -1.8),
      pval_synth = c(0.01, 0.03),
      fdr_synth = c(0.05, 0.1)
    ),
    exp2 = data.frame(
      GENE = c("Trp53", "Brca1"),
      numObs = c(4, 5),
      normZ = c(1.5, -2.1),
      pval_synth = c(0.02, 0.01),
      fdr_synth = c(0.08, 0.05)
    )
  )

  # With pattern = NULL, should use all items and full names as prefixes
  result <- BuildDrugZdz(drugz_list, pattern = NULL)

  expect_true(any(grepl("exp1", names(result))))
  expect_true(any(grepl("exp2", names(result))))
})
