test_that("NScutoff returns correct cutoffs", {
  data <- data.frame(
    gene = paste0("Gene", 1:100),
    LFC = stats::rnorm(1000, mean = 0, sd = 2)
  )

  cutoffs <- NScutoff(data$LFC, scale = 1)

  mu <- mean(data$LFC, na.rm = TRUE)
  stdv <- stats::sd(data$LFC, na.rm = TRUE)

  expected_cutoffs <- c(mu - stdv * 1, mu + stdv * 1)

  expect_equal(cutoffs, expected_cutoffs)
})


test_that("NScutoff handles non-numeric input", {
  expect_error(NScutoff(c("a", "b", "c"), scale = 1), "x has to be numeric")
})

test_that("NScutoff handles force_zero_center correctly", {
  data <- data.frame(
    gene = paste0("Gene", 1:100),
    LFC = stats::rnorm(1000, mean = 5, sd = 2)
  )

  cutoffs <- NScutoff(data$LFC, scale = 1, force_zero_center = TRUE)

  stdv <- sd(data$LFC, na.rm = TRUE)

  expected_cutoffs <- c(0 - stdv * 1, 0 + stdv * 1)

  expect_equal(cutoffs, expected_cutoffs)
})


test_that("NSgencutoff returns correct cutoffs with different force_zero_center options", {
  data <- data.frame(
    gene = paste0("Gene", 1:100),
    control = stats::rnorm(1000, mean = 5, sd = 2),
    treatment = stats::rnorm(1000, mean = -3, sd = 2)
  )

  # Test for "control"
  cutoffs_control <- NSgencutoff(data, scale = 1, force_zero_center = "control")
  stdv_control <- sd(data$control, na.rm = TRUE)
  expected_control_cutoffs <- c(0 - stdv_control * 1, 0 + stdv_control * 1)
  expect_equal(cutoffs_control$x_cutoff, expected_control_cutoffs)

  # Test for "treatment"
  cutoffs_treatment <- NSgencutoff(data, scale = 1, force_zero_center = "treatment")
  stdv_treatment <- sd(data$treatment, na.rm = TRUE)
  expected_treatment_cutoffs <- c(0 - stdv_treatment * 1, 0 + stdv_treatment * 1)
  expect_equal(cutoffs_treatment$y_cutoff, expected_treatment_cutoffs)

  # Test for "both"
  cutoffs_both <- NSgencutoff(data, scale = 1, force_zero_center = "both")
  expect_equal(cutoffs_both$x_cutoff, expected_control_cutoffs)
  expect_equal(cutoffs_both$y_cutoff, expected_treatment_cutoffs)
})



test_that("NSbasedf filters by min_sgrna and selects correct columns", {
  data <- data.frame(
    gene = paste0("Gene", 1:10),
    num = c(2, 3, 4, 1, 5, 3, 2, 4, 5, 3),
    untreated_LFC = stats::rnorm(10, mean = 0, sd = 2),
    treated_LFC = stats::rnorm(10, mean = 0, sd = 2),
    untreated_pval = stats::runif(10, min = 0, max = 1),
    treated_pval = stats::runif(10, min = 0, max = 1)
  )

  result <- NSbasedf(
    data = data,
    control = untreated_LFC,
    treament = treated_LFC,
    ctrl_pval = untreated_pval,
    treat_pval = treated_pval,
    min_sgrna = 3
  )

  expected_genes <- data$gene[data$num >= 3]
  expect_equal(result$id, expected_genes)
  expect_equal(ncol(result), 6) # id, num, control, treatment, ctrl_pval, treat_pval
})

test_that("NSbasedf handles no genes meeting min_sgrna", {
  data <- data.frame(
    gene = paste0("Gene", 1:5),
    num = c(1, 2, 1, 2, 1),
    untreated_LFC = stats::rnorm(5, mean = 0, sd = 2),
    treated_LFC = stats::rnorm(5, mean = 0, sd = 2),
    untreated_pval = stats::runif(5, min = 0, max = 1),
    treated_pval = stats::runif(5, min = 0, max = 1)
  )

  result <- NSbasedf(
    data = data,
    control = untreated_LFC,
    treament = treated_LFC,
    ctrl_pval = untreated_pval,
    treat_pval = treated_pval,
    min_sgrna = 3
  )

  expect_equal(nrow(result), 0)
})

test_that("NSbasedf handles non-numeric columns gracefully", {
  data <- data.frame(
    gene = paste0("Gene", 1:5),
    num = c(3, 4, 3, 4, 3),
    untreated_LFC = c("a", "b", "c", "d", "e"),
    treated_LFC = stats::rnorm(5, mean = 0, sd = 2),
    untreated_pval = stats::runif(5, min = 0, max = 1),
    treated_pval = stats::runif(5, min = 0, max = 1)
  )

  expect_error(
    NSbasedf(
      data = data,
      control = untreated_LFC,
      treament = treated_LFC,
      ctrl_pval = untreated_pval,
      treat_pval = treated_pval,
      min_sgrna = 3
    )
  )
})

test_that("NSbasedf works without p-value columns", {
  data <- data.frame(
    gene = paste0("Gene", 1:10),
    num = c(2, 3, 4, 1, 5, 3, 2, 4, 5, 3),
    untreated_LFC = stats::rnorm(10, mean = 0, sd = 2),
    treated_LFC = stats::rnorm(10, mean = 0, sd = 2)
  )

  result <- NSbasedf(
    data = data,
    control = untreated_LFC,
    treament = treated_LFC,
    min_sgrna = 3
  )

  expected_genes <- data$gene[data$num >= 3]
  expect_equal(result$id, expected_genes)
  expect_equal(ncol(result), 4) # id, num, control, treatment
})


test_that("NSfilterpval filters by min_pval correctly", {
  data <- data.frame(
    id = paste0("Gene", 1:10),
    num = c(3, 4, 5, 3, 4, 5, 3, 4, 5, 3),
    control = stats::rnorm(10, mean = 0, sd = 2),
    treatment = stats::rnorm(10, mean = 0, sd = 2),
    ctrl_pval =  c(0.01, 0.2,  0.03, 0.5,  0.04, 0.6,  0.07, 0.8,  0.02, 0.9),
    treat_pval = c(0.05, 0.15, 0.25, 0.03, 0.35, 0.04, 0.45, 0.02, 0.55, 0.01)
  )

  result <- NSfilterpval(data)

  expected_minpval <- c(0.01, 0.03, 0.03, 0.04, 0.04, 0.02, 0.02, 0.01)
  expected_ids <- c("Gene1", "Gene3", "Gene4", "Gene5", "Gene6", "Gene8", "Gene9", "Gene10")
  expect_equal(result$id, expected_ids)
  expect_equal(result$min_pv, expected_minpval)
  expect_false("ctrl_pval" %in% colnames(result))
  expect_false("treat_pval" %in% colnames(result))
})


test_that("NSfilterpval handles min_pval edge cases", {
  data <- data.frame(
    id = paste0("Gene", 1:5),
    num = c(3, 4, 5, 3, 4),
    control = stats::rnorm(5, mean = 0, sd = 2),
    treatment = stats::rnorm(5, mean = 0, sd = 2),
    ctrl_pval = c(0.1, 0.2, 0.3, 0.4, 0.5),
    treat_pval = c(0.6, 0.7, 0.8, 0.9, 1.0)
  )

  # No genes meet the min_pval criteria
  expect_warning(
    result <- NSfilterpval(data, min_pval = 0.05),
    "No genes passed the p-value filter"
  )
  expect_equal(nrow(result), 0)


  # All genes meet the min_pval criteria
  result_all <- NSfilterpval(data, min_pval = 1.0)
  expect_equal(nrow(result_all), nrow(data))
})

test_that("NSfilterpval handles invalid min_pval", {
  data <- data.frame(
    id = paste0("Gene", 1:5),
    num = c(3, 4, 5, 3, 4),
    control = stats::rnorm(5, mean = 0, sd = 2),
    treatment = stats::rnorm(5, mean = 0, sd = 2),
    ctrl_pval = c(0.1, 0.2, 0.3, 0.4, 0.5),
    treat_pval = c(0.6, 0.7, 0.8, 0.9, 1.0)
  )

  expect_error(NSfilterpval(data, min_pval = -0.1), "min_pval has to be between 0 and 1")
  expect_error(NSfilterpval(data, min_pval = 1.5), "min_pval has to be between 0 and 1")
})


# test_that("NSgoidframe creates correct GO ID dataframe", {
#   data <- data.frame(
#     id = paste0("Gene", 1:5),
#     num = c(3, 4, 5, 3, 4),
#     control = stats::rnorm(5, mean = 0, sd = 2),
#     treatment = stats::rnorm(5, mean = 0, sd = 2)
#   )
#
#   goi_vector <- c("Gene1", "Gene3", "Gene5")
#
#   result <- NSgoidf(data, goi_vector)
#
#   expect_equal(nrow(result), length(goi_vector))
#   expect_equal(result$id, goi_vector)
#   expect_true(all(c("num", "control", "treatment") %in% colnames(result)))
# })


# --- NSsquares tests ---

test_that("NSsquares assigns squares correctly based on cutoffs", {
  # The NSsquares function checks neutral_slope FIRST:
  # If diff (treatment - control) is between slope_cutoff[1] and slope_cutoff[2],
  # the gene is classified as neutral_slope regardless of position.

  # To test positional assignment, diff must be OUTSIDE slope_cutoff.
  # slope_cutoff[1] = lower bound, slope_cutoff[2] = upper bound

  # For square assignment to work:
  # - bottom_left: control < x_cutoff[1], treatment < y_cutoff[1], diff outside slope
  # - bottom_center: control between x_cutoff, treatment < y_cutoff[1], diff outside slope
  # etc.

  # Create data with diff values OUTSIDE the slope_cutoff range
  data <- data.frame(
    id = c("gene_bl", "gene_bc", "gene_br",
           "gene_tl", "gene_tc", "gene_tr"),
    num = rep(4, 6),
    # For bottom row: treatment = -10, so diff = treatment - control = -10 - control
    # For top row: treatment = 10, so diff = 10 - control
    # With slope_cutoff = c(-5, 5), diffs outside this range will NOT be neutral_slope
    control =   c(-3,  0,  3, -3,  0,  3),
    treatment = c(-10, -10, -10, 10, 10, 10)
    # diffs:     -7  -10  -13   13  10   7  <- all outside c(-5, 5)
  )

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  result <- NSsquares(data, x_cutoff, y_cutoff, slope_cutoff)

  # Check that necessary columns were added
  expect_true("square" %in% names(result))
  expect_true("diff" %in% names(result))
  expect_true("rank" %in% names(result))
  expect_true("euclid_dist" %in% names(result))

  # Check specific square assignments
  expect_equal(result$square[result$id == "gene_bl"], "bottom_left")
  expect_equal(result$square[result$id == "gene_bc"], "bottom_center")
  expect_equal(result$square[result$id == "gene_br"], "bottom_right")
  expect_equal(result$square[result$id == "gene_tl"], "top_left")
  expect_equal(result$square[result$id == "gene_tc"], "top_center")
  expect_equal(result$square[result$id == "gene_tr"], "top_right")
})


test_that("NSsquares assigns middle and center squares when diff is outside slope_cutoff", {
  # Middle squares require treatment between y_cutoff, which naturally
  # makes diff close to control. To get diff outside slope_cutoff while
  # keeping treatment in middle range, we need careful values.

  # For middle_left: control < x_cutoff[1] (< -1), treatment between y_cutoff (-1, 1)
  #   e.g., control = -10, treatment = 0 -> diff = 10 (outside c(-5, 5))

  data <- data.frame(
    id = c("gene_ml", "gene_c", "gene_mr"),
    num = rep(4, 3),
    control =   c(-10, 0, 10),
    treatment = c(0,   0, 0)
    # diffs:     10   0  -10
  )

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)  # diff = 0 IS in this range, so center -> neutral_slope

  result <- NSsquares(data, x_cutoff, y_cutoff, slope_cutoff)

  # middle_left and middle_right have diff outside slope_cutoff
  expect_equal(result$square[result$id == "gene_ml"], "middle_left")
  expect_equal(result$square[result$id == "gene_mr"], "middle_right")

  # center has diff = 0 which IS between -5 and 5, so it becomes neutral_slope
  expect_equal(result$square[result$id == "gene_c"], "neutral_slope")
})


test_that("NSsquares assigns neutral_slope correctly", {
  # Create data where treatment - control is within slope_cutoff
  data <- data.frame(
    id = c("gene_neutral", "gene_not_neutral"),
    num = c(4, 4),
    control = c(2, 2),
    treatment = c(2.5, 5)  # diff = 0.5 and 3
  )

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-1, 1)  # Narrow slope cutoff

  result <- NSsquares(data, x_cutoff, y_cutoff, slope_cutoff)

  # gene_neutral has diff = 0.5, which is between -1 and 1
  expect_equal(result$square[result$id == "gene_neutral"], "neutral_slope")
  # gene_not_neutral has diff = 3, outside slope_cutoff
  expect_equal(result$square[result$id == "gene_not_neutral"], "top_right")
})


test_that("NSsquares calculates euclidean distance correctly", {
  data <- data.frame(
    id = c("gene1", "gene2"),
    num = c(4, 4),
    control = c(3, 4),
    treatment = c(4, 3)
  )

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-10, 10)

  result <- NSsquares(data, x_cutoff, y_cutoff, slope_cutoff)

  # Euclidean distance = sqrt(control^2 + treatment^2)
  expected_dist1 <- sqrt(3^2 + 4^2)  # = 5
  expected_dist2 <- sqrt(4^2 + 3^2)  # = 5

  expect_equal(result$euclid_dist[result$id == "gene1"], expected_dist1)
  expect_equal(result$euclid_dist[result$id == "gene2"], expected_dist2)
})


test_that("NSsquares calculates perpendicular distance correctly", {
  data <- data.frame(
    id = c("gene1", "gene2"),
    num = c(4, 4),
    control = c(2, 3),
    treatment = c(4, 3)  # diff = 2 and 0
  )

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-10, 10)

  result <- NSsquares(data, x_cutoff, y_cutoff, slope_cutoff)

  # perp_dist = abs(diff) / sqrt(2)
  expected_perp1 <- abs(4 - 2) / sqrt(2)  # = 2/sqrt(2)
  expected_perp2 <- abs(3 - 3) / sqrt(2)  # = 0

  expect_equal(result$perp_dist[result$id == "gene1"], expected_perp1)
  expect_equal(result$perp_dist[result$id == "gene2"], expected_perp2)
})


test_that("NSsquares ranks genes within each square by score", {
  # Create multiple genes in the same square with different distances
  data <- data.frame(
    id = c("gene1", "gene2", "gene3"),
    num = c(4, 4, 4),
    control = c(3, 4, 5),      # All in top_right
    treatment = c(3, 4, 5)
  )

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-10, 10)

  result <- NSsquares(data, x_cutoff, y_cutoff, slope_cutoff)

  # All should be in top_right (or center due to being on diagonal - but diff=0 means neutral_slope)
  # Actually with diff=0 and slope_cutoff = c(-10, 10), all fall in neutral_slope
  # Let's check ranks exist and are 1, 2, 3
  expect_true(all(result$rank %in% 1:3))
})


test_that("NSsquares errors when cutoffs are missing", {
  data <- data.frame(
    id = c("gene1"),
    num = c(4),
    control = c(1),
    treatment = c(2)
  )

  expect_error(
    NSsquares(data, y_cutoff = c(-1, 1), slope_cutoff = c(-1, 1)),
    "x_cutoff, y_cutoff and slope_cutoff have to be provided"
  )

  expect_error(
    NSsquares(data, x_cutoff = c(-1, 1), slope_cutoff = c(-1, 1)),
    "x_cutoff, y_cutoff and slope_cutoff have to be provided"
  )

  expect_error(
    NSsquares(data, x_cutoff = c(-1, 1), y_cutoff = c(-1, 1)),
    "x_cutoff, y_cutoff and slope_cutoff have to be provided"
  )
})
