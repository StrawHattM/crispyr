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
