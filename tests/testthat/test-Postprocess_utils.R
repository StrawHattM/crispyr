# Tests for post-processing utility functions: extractNSdata and where_is

# Helper function to create a mock Nine Squares graph object
create_mock_ns_graph <- function(genes = paste0("Gene", 1:10)) {
  data <- data.frame(
    id = genes,
    num = rep(4, length(genes)),
    control = c(-3, -2, 0, 1, 2, 3, -2, 0, 2, 0),
    treatment = c(-2, 3, -3, 0, 3, 2, 0, 0, -2, 1),
    square = c("bottom_left", "top_left", "bottom_center", "center",
               "top_right", "top_right", "middle_left", "center",
               "bottom_right", "center"),
    rank = c(1, 1, 1, 1, 1, 2, 1, 2, 1, 3),
    diff = c(1, 5, -3, -1, 1, -1, 2, 0, -4, 1),
    euclid_dist = sqrt(c(13, 13, 9, 1, 13, 13, 4, 0, 8, 1)),
    perp_dist = abs(c(1, 5, -3, -1, 1, -1, 2, 0, -4, 1)) / sqrt(2),
    score = 1:10
  )

  ggplot2::ggplot(data, ggplot2::aes(x = control, y = treatment)) +
    ggplot2::geom_point()
}


# --- extractNSdata tests ---

test_that("extractNSdata extracts full data from single graph", {
  graph <- create_mock_ns_graph()

  result <- extractNSdata(graph, which = "data")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 10)
  expect_true("id" %in% names(result))
  expect_true("square" %in% names(result))
})


test_that("extractNSdata extracts genes from specific squares", {
  graph <- create_mock_ns_graph()

  # Extract genes from top_right
  result <- extractNSdata(graph, which = "top_right")

  expect_type(result, "character")
  expect_true("Gene5" %in% result)
  expect_true("Gene6" %in% result)
  expect_equal(length(result), 2)
})


test_that("extractNSdata extracts genes from multiple squares", {
  graph <- create_mock_ns_graph()

  result <- extractNSdata(graph, which = c("top_right", "bottom_left"))

  expect_type(result, "character")
  expect_true("Gene5" %in% result)  # top_right
  expect_true("Gene1" %in% result)  # bottom_left
})


test_that("extractNSdata errors when 'data' combined with other options", {
  graph <- create_mock_ns_graph()

  expect_error(
    extractNSdata(graph, which = c("data", "top_right")),
    "If 'data' is selected, it must be the only option selected"
  )
})


test_that("extractNSdata handles list of graphs", {
  graph1 <- create_mock_ns_graph(paste0("GeneA", 1:5))
  graph2 <- create_mock_ns_graph(paste0("GeneB", 1:5))

  graph_list <- list(
    comparison1 = graph1,
    comparison2 = graph2
  )

  result <- extractNSdata(graph_list, which = "data")

  expect_type(result, "list")
  expect_equal(length(result), 2)
  expect_true(all(c("comparison1", "comparison2") %in% names(result)))
  expect_s3_class(result$comparison1, "data.frame")
})


test_that("extractNSdata extracts squares from list of graphs", {
  graph1 <- create_mock_ns_graph()
  graph2 <- create_mock_ns_graph()

  graph_list <- list(
    screen1 = graph1,
    screen2 = graph2
  )

  result <- extractNSdata(graph_list, which = "top_right")

  expect_type(result, "list")
  expect_true(all(sapply(result, is.character)))
})


test_that("extractNSdata validates 'which' parameter", {
  graph <- create_mock_ns_graph()

  expect_error(
    extractNSdata(graph, which = "invalid_square"),
    "'arg' should be one of"
  )
})


# --- where_is tests ---

test_that("where_is finds single gene in dataframe", {
  data <- data.frame(
    id = paste0("Gene", 1:5),
    square = c("top_left", "center", "bottom_right", "top_center", "middle_left")
  )

  result <- where_is("Gene1", data)

  expect_equal(result, "top_left")
})


test_that("where_is finds multiple genes in dataframe", {
  data <- data.frame(
    id = paste0("Gene", 1:5),
    square = c("top_left", "center", "bottom_right", "top_center", "middle_left")
  )

  result <- where_is(c("Gene1", "Gene3", "Gene5"), data)

  expect_type(result, "character")
  expect_named(result, c("Gene1", "Gene3", "Gene5"))
  expect_equal(result["Gene1"], c(Gene1 = "top_left"))
  expect_equal(result["Gene3"], c(Gene3 = "bottom_right"))
})


test_that("where_is returns NA for genes not found", {
  data <- data.frame(
    id = paste0("Gene", 1:3),
    square = c("top_left", "center", "bottom_right")
  )

  result <- where_is(c("Gene1", "Gene99"), data)

  expect_equal(result["Gene1"], c(Gene1 = "top_left"))
  expect_true(is.na(result["Gene99"]))
})


test_that("where_is extracts data from ggplot object", {
  graph <- create_mock_ns_graph()

  result <- where_is("Gene1", graph)

  expect_equal(result, "bottom_left")
})


test_that("where_is handles list of ggplot objects", {
  graph1 <- create_mock_ns_graph()
  graph2 <- create_mock_ns_graph()

  # Modify graph2 data to have Gene1 in a different square
  graph2$data$square[1] <- "top_center"

  graph_list <- list(
    screen1 = graph1,
    screen2 = graph2
  )

  result <- where_is("Gene1", graph_list)

  expect_type(result, "list")
  expect_equal(names(result), c("screen1", "screen2"))
  expect_equal(result$screen1, "bottom_left")
  expect_equal(result$screen2, "top_center")
})


test_that("where_is handles list of dataframes", {
  df1 <- data.frame(
    id = paste0("Gene", 1:3),
    square = c("top_left", "center", "bottom_right")
  )

  df2 <- data.frame(
    id = paste0("Gene", 1:3),
    square = c("bottom_left", "top_right", "center")
  )

  df_list <- list(
    condition1 = df1,
    condition2 = df2
  )

  result <- where_is("Gene1", df_list)

  expect_type(result, "list")
  expect_equal(result$condition1, "top_left")
  expect_equal(result$condition2, "bottom_left")
})


test_that("where_is handles named list of character vectors", {
  gene_lists <- list(
    top_genes = c("Gene1", "Gene2", "Gene3"),
    bottom_genes = c("Gene4", "Gene5"),
    other = c("Gene6", "Gene7", "Gene8")
  )

  result <- where_is("Gene2", gene_lists)

  expect_equal(result, "top_genes")
})


test_that("where_is returns character(0) for gene not in any list", {
  gene_lists <- list(
    group1 = c("Gene1", "Gene2"),
    group2 = c("Gene3", "Gene4")
  )

  result <- where_is("Gene99", gene_lists)

  expect_equal(result, character(0))
})


test_that("where_is return_full parameter works with list input", {
  gene_lists <- list(
    group1 = c("Gene1", "Gene2", "Gene3"),
    group2 = c("Gene4", "Gene5")
  )

  # Without return_full
  result_names <- where_is("Gene2", gene_lists, return_full = FALSE)
  expect_equal(result_names, "group1")

  # With return_full
  result_full <- where_is("Gene2", gene_lists, return_full = TRUE)
  expect_type(result_full, "list")
  expect_equal(result_full$group1, c("Gene1", "Gene2", "Gene3"))
})


test_that("where_is handles multiple genes with list of graphs", {
  graph1 <- create_mock_ns_graph()
  graph2 <- create_mock_ns_graph()

  graph_list <- list(
    exp1 = graph1,
    exp2 = graph2
  )

  result <- where_is(c("Gene1", "Gene5"), graph_list)

  expect_type(result, "list")
  expect_equal(names(result), c("exp1", "exp2"))

  # Each element should be a named character vector
  expect_named(result$exp1, c("Gene1", "Gene5"))
})