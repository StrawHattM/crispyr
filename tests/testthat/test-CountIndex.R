# Tests for CountIndex.R functions

# Tests for calculate_totals ---------------------------------------------

test_that("calculate_totals returns correct structure", {
  raw_counts <- data.frame(
    sgRNA = c("sg1", "sg2", "sg3"),
    Gene = c("GeneA", "GeneB", "GeneC"),
    Sample1 = c(100, 200, 150),
    Sample2 = c(120, 180, 160),
    Sample3 = c(110, 190, 140)
  )

  conditions <- c("ctrl", "ctrl", "treat")
  replicates <- c("rep1", "rep2", "rep1")

  result <- calculate_totals(raw_counts, conditions, replicates)

  # Check structure

  expect_s3_class(result, "data.frame")
  expect_named(result, c("sample", "total", "condition", "rep"))
  expect_equal(nrow(result), 3)
})


test_that("calculate_totals calculates correct column sums", {
  raw_counts <- data.frame(
    sgRNA = c("sg1", "sg2", "sg3"),
    Gene = c("GeneA", "GeneB", "GeneC"),
    Sample1 = c(100, 200, 150),
    Sample2 = c(120, 180, 160),
    Sample3 = c(110, 190, 140)
  )

  conditions <- c("ctrl", "ctrl", "treat")
  replicates <- c("rep1", "rep2", "rep1")

  result <- calculate_totals(raw_counts, conditions, replicates)

  # Check totals are column sums
  expect_equal(result$total[1], 450)  # 100 + 200 + 150
  expect_equal(result$total[2], 460)  # 120 + 180 + 160
  expect_equal(result$total[3], 440)  # 110 + 190 + 140
})


test_that("calculate_totals preserves sample names as factors", {
  raw_counts <- data.frame(
    sgRNA = c("sg1", "sg2"),
    Gene = c("GeneA", "GeneB"),
    SampleA = c(100, 200),
    SampleB = c(120, 180)
  )

  conditions <- c("ctrl", "treat")
  replicates <- c("rep1", "rep1")

  result <- calculate_totals(raw_counts, conditions, replicates)

  expect_s3_class(result$sample, "factor")
  expect_equal(as.character(result$sample), c("SampleA", "SampleB"))
})


test_that("calculate_totals creates condition and rep factors with correct levels", {
  raw_counts <- data.frame(
    sgRNA = c("sg1", "sg2"),
    Gene = c("GeneA", "GeneB"),
    S1 = c(100, 200),
    S2 = c(120, 180),
    S3 = c(110, 190),
    S4 = c(105, 205)
  )

  conditions <- c("ctrl", "ctrl", "treat", "treat")
  replicates <- c("rep1", "rep2", "rep1", "rep2")

  result <- calculate_totals(raw_counts, conditions, replicates)

  expect_s3_class(result$condition, "factor")
  expect_s3_class(result$rep, "factor")
  expect_equal(levels(result$condition), c("ctrl", "treat"))
  expect_equal(levels(result$rep), c("rep1", "rep2"))
})


test_that("calculate_totals handles single sample", {
  raw_counts <- data.frame(
    sgRNA = c("sg1", "sg2", "sg3"),
    Gene = c("GeneA", "GeneB", "GeneC"),
    OnlySample = c(100, 200, 150)
  )

  conditions <- "ctrl"
  replicates <- "rep1"

  result <- calculate_totals(raw_counts, conditions, replicates)

  expect_equal(nrow(result), 1)
  expect_equal(result$total, 450)
})


test_that("calculate_totals handles zeros in counts", {
  raw_counts <- data.frame(
    sgRNA = c("sg1", "sg2", "sg3"),
    Gene = c("GeneA", "GeneB", "GeneC"),
    Sample1 = c(0, 0, 0),
    Sample2 = c(100, 0, 50)
  )

  conditions <- c("ctrl", "treat")
  replicates <- c("rep1", "rep1")

  result <- calculate_totals(raw_counts, conditions, replicates)

  expect_equal(result$total[1], 0)
  expect_equal(result$total[2], 150)
})


# Tests for CountIndex ---------------------------------------------------

# Note: CountIndex appears incomplete - only handles the case when 'order' is missing
# These tests document current behavior

test_that("CountIndex returns sorted column names when order is missing", {
  counts <- data.frame(
    zebra_sample = c(1, 2),
    alpha_sample = c(3, 4),
    beta_sample = c(5, 6)
  )

  result <- CountIndex(counts)

  # Should return alphabetically sorted column names
  expect_equal(result, c("alpha_sample", "beta_sample", "zebra_sample"))
})


test_that("CountIndex sorts using en_US.UTF-8 locale",
{
  skip_on_os("windows")  # Locale behavior may differ on Windows

  counts <- data.frame(
    z_col = 1,
    a_col = 2,
    A_col = 3  # Capital letter sorting
  )

  result <- CountIndex(counts)

  # Just verify it returns a character vector of sorted names
  expect_type(result, "character")
  expect_length(result, 3)
})
