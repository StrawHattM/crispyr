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
