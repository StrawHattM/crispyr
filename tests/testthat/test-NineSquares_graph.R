# Tests for Nine Squares graph-building functions

# Helper function to create test data with square assignments
create_test_data_with_squares <- function() {
  data <- data.frame(
    id = paste0("Gene", 1:20),
    num = rep(4, 20),
    control = c(-3, -3, -3, -3, 0, 0, 0, 0, 3, 3, 3, 3, -2, 2, -1, 1, 0, 0.5, -0.5, 0),
    treatment = c(-3, 0, 3, -2, -3, 0, 3, -2, -3, 0, 3, 2, 3, -3, 0, 0, 1, -1, 2, -2)
  )

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  NSsquares(data, x_cutoff, y_cutoff, slope_cutoff)
}


# --- NSbasegraph tests ---

test_that("NSbasegraph returns a ggplot object", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff)

  expect_s3_class(graph, "gg")
  expect_s3_class(graph, "ggplot")
})


test_that("NSbasegraph includes correct layers", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff)

  # Check that data is attached
  expect_equal(nrow(graph$data), nrow(data))

  # Check that the graph has layers (points, vlines, hlines, ablines)
  expect_gt(length(graph$layers), 0)
})


test_that("NSbasegraph errors on invalid legend parameter", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  expect_error(
    NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff, legend = "yes"),
    "legend needs to be TRUE or FALSE"
  )
})


test_that("NSbasegraph handles size_var parameter", {
  data <- create_test_data_with_squares()
  data$my_size <- runif(nrow(data), 1, 5)

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  # With size_var - should not error
  graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff, size_var = "my_size")

  expect_s3_class(graph, "ggplot")
})


test_that("NSbasegraph respects alpha and shape parameters", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff,
                       alpha = 0.8, shape = 16, size = 3)

  expect_s3_class(graph, "ggplot")
})


# --- NSaxislabels tests ---

test_that("NSaxislabels adds axis labels to graph", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  base_graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff)

  labeled_graph <- NSaxislabels(base_graph, xlab = "Control LFC", ylab = "Treatment LFC")

  expect_s3_class(labeled_graph, "ggplot")
  expect_equal(labeled_graph$labels$x, "Control LFC")
  expect_equal(labeled_graph$labels$y, "Treatment LFC")
})


# --- NStitle tests ---

test_that("NStitle adds title to graph", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  base_graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff)

  titled_graph <- NStitle(base_graph, title = "My Nine Squares Plot")

  expect_s3_class(titled_graph, "ggplot")
  expect_equal(titled_graph$labels$title, "My Nine Squares Plot")
})


# --- NSaddtoplabels tests ---

test_that("NSaddtoplabels adds text labels to graph", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  base_graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff)

  labeled_graph <- NSaddtoplabels(
    base_graph, data,
    groups_labeled = c("top_center", "bottom_center"),
    top_labeled = 3
  )

  expect_s3_class(labeled_graph, "ggplot")
  # Check that a text layer was added
  layer_count_before <- length(base_graph$layers)
  layer_count_after <- length(labeled_graph$layers)

  expect_gt(layer_count_after, layer_count_before)
})


test_that("NSaddtoplabels handles different groups_labeled", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  base_graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff)

  # Test with corner squares
  labeled_graph <- NSaddtoplabels(
    base_graph, data,
    groups_labeled = c("top_left", "bottom_right", "top_right", "bottom_left"),
    top_labeled = 2
  )

  expect_s3_class(labeled_graph, "ggplot")
})


# --- NSaddgoi tests ---

test_that("NSaddgoi highlights genes of interest", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  base_graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff)

  goi_graph <- NSaddgoi(
    data = data,
    graph = base_graph,
    goi_list = c("Gene1", "Gene5", "Gene10")
  )

  expect_s3_class(goi_graph, "ggplot")

  # Check that layers were added
  expect_gt(length(goi_graph$layers), length(base_graph$layers))
})


test_that("NSaddgoi works with goi_auto = TRUE", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  base_graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff)

  goi_graph <- NSaddgoi(
    data = data,
    graph = base_graph,
    goi_list = c("Gene1", "Gene5"),
    goi_auto = TRUE
  )

  expect_s3_class(goi_graph, "ggplot")
})


test_that("NSaddgoi respects goi_label_type parameter", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  base_graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff)

  # Test with "text" label type
  goi_graph_text <- NSaddgoi(
    data = data,
    graph = base_graph,
    goi_list = c("Gene1"),
    goi_label_type = "text"
  )

  expect_s3_class(goi_graph_text, "ggplot")

  # Test with "label" label type (default)
  goi_graph_label <- NSaddgoi(
    data = data,
    graph = base_graph,
    goi_list = c("Gene1"),
    goi_label_type = "label"
  )

  expect_s3_class(goi_graph_label, "ggplot")
})


test_that("NSaddgoi handles custom color and size parameters", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  base_graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff)

  goi_graph <- NSaddgoi(
    data = data,
    graph = base_graph,
    goi_list = c("Gene1"),
    goi_color = "red",
    goi_fill = "blue",
    goi_size = 4,
    goi_shape = 23,
    goi_label_color = "darkred",
    goi_label_size = 5
  )

  expect_s3_class(goi_graph, "ggplot")
})


test_that("NSaddgoi can extract data from graph if data not provided", {
  data <- create_test_data_with_squares()

  x_cutoff <- c(-1, 1)
  y_cutoff <- c(-1, 1)
  slope_cutoff <- c(-5, 5)

  base_graph <- NSbasegraph(data, x_cutoff, y_cutoff, slope_cutoff)

  # Don't provide data argument - should extract from graph
  goi_graph <- NSaddgoi(
    graph = base_graph,
    goi_list = c("Gene1", "Gene5")
  )

  expect_s3_class(goi_graph, "ggplot")
})