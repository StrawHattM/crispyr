
# System setup: adding Rprofile data --------------------------------------

## I added this code to Rprofile


# if (interactive()) {
#   suppressMessages(require(devtools))
# }
#
# options(
#   usethis.description = list(
#     "Authors@R" = utils::person(
#       "Martín", "González Fernández",
#       email = "martin.gonzalezfernandez@unibe.ch",
#       role = c("aut", "cre"),
#       comment = c(ORCID = "0009-0006-8308-5493")
#     ),
#     License = "MIT + file LICENSE"
#   )
# )



# UseThese -----------------------------------------------------------------

use_github("StrawHattM/crispyr")

use_mit_license()

use_testthat()

usethis::use_package("dplyr", min_version = TRUE)
usethis::use_package("ggplot2", min_version = TRUE)
usethis::use_package("ggrepel", min_version = TRUE)
usethis::use_package("ggnewscale", min_version = TRUE)
usethis::use_package("readr", min_version = TRUE)
usethis::use_package("stringr", min_version = TRUE)
usethis::use_package("tidyr", min_version = TRUE)





usethis::use_import_from("magrittr", "%>%")
usethis::use_import_from("rlang", ".data")
usethis::use_import_from("purrr", "reduce")
usethis::use_import_from("withr", "with_locale")


# First function: %nin% ---------------------------------------------------

usethis::use_r("nin")

## here I added %nin% code to the R file and saved it.

load_all() # should make the %nin% function available

## we can check it works after loading
c("A", "B", "C", "D", "E", "F") %nin% c("B", "A", "D", "E")

exists("%nin%", where = globalenv(), inherits = FALSE)
## has to throw FALSE because we're in the development environment, not global

check() # global check for getting used to it

## When in the function, create a roxygen2 skeleton (Ctrl+Alt+Shift+R) to document
## the function and save it. After, call document()
document()

## This is weird because %nin% shouldn't be a function, but a logical operator
?`%nin%` # necessary to add `` to escape the special characters.

install()

## Make tests, again weird because of the operator nature

usethis::use_test("grapes-nin-grapes")

test()


# NineSquares helpers for the big one ------------------------------------------------

usethis::use_r("NineSquares_helpers")

document()
load_all()



## This goes into the examples of the function, but it's good to check that it works tbh
data <- data.frame(gene = paste0("Gene", 1:100),
                   LFC = rnorm(1000, mean = 5, sd = 2))
cutoffs <- NScutoff(data$LFC, scale = 1)



usethis::use_test("NineSquares_helpers")
test()

# checking behaviour of NSbasedf

set.seed(42)

data <- data.frame(
  gene = paste0("Gene", 1:20000),
  num = sample(1:10, 20000, replace = TRUE),
  untreated_LFC = stats::rnorm(20000, mean = 0, sd = 3),
  treated_LFC = stats::rnorm(20000, mean = 0, sd = 3),
  untreated_pval = stats::runif(20000, min = 0, max = 1),
  treated_pval = stats::runif(20000, min = 0, max = 1)
)

base_data <-
  NSbasedf(data, control = untreated_LFC,
           treament = treated_LFC,
           ctrl_pval = untreated_pval,
           treat_pval = treated_pval,
           min_sgrna = 2)

cutoffs <-
  NSgencutoff(base_data, scale = 2, force_zero_center = "none")

base_data <-
  NSsquares(base_data,
            x_cutoff = cutoffs$x_cutoff,
            y_cutoff = cutoffs$y_cutoff,
            slope_cutoff = cutoffs$slope_cutoff)

base_graph<-
  NSbasegraph(data = base_data,
              alpha = 0.4,
              shape = 21,
              size = 2,
              x_cutoff = cutoffs$x_cutoff,
              y_cutoff = cutoffs$y_cutoff,
              slope_cutoff = cutoffs$slope_cutoff)
# graph_test

# base_graph_toplabels <-
  NSaddtoplabels(base_graph, data = base_data, top_labeled = 5)


graph_test %>%
  NSaddgoi(data = base_data,
           graph = graph_test,
           goi_list = c("Gene2994", "Gene405", "Gene1450", "Gene592", "Gene14261"),
           goi_auto = TRUE)


# NineSquares assembly ----------------------------------------------------


usethis::use_r("NineSquares")

NineSquares(data,
            untreated_LFC,
            treated_LFC,
            untreated_pval,
            treated_pval,
            # min_pval = 0.25,
            scale = 2,
            shape = 15,
            top_labeled = 10,
            xlab = "Stupidly Awesome",
            ylab = "Awesomely Stupid",
            title = "No way you're doing this?",
            legend = TRUE,
            filename = "testing_whynot",
            groups_labeled = c("top_center", "bottom_center"),
            goi = c("Gene2994", "Gene405", "Gene1450", "Gene592", "Gene14261"))


## Works!


# Parsing raw counts for groups -------------------------------------------

## Since everything is so very dependent on number of columns and column grouping,
## there should be some starting function to establish groupings that then can be
## inherited by other functions, a metadata of sorts. Could be a SummarizedExperiment.


# Quality Control ---------------------------------------------------------


## Count index

usethis::use_r("CountIndex")










# Import and process RRA --------------------------------------------------

usethis::use_r("ImportRRA")

# Import and process MLE --------------------------------------------------

usethis::use_r("ImportMLE")
