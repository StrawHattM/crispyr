
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

use_mit_licuseense()

use_testthat()

usethis::use_package("dplyr", min_version = TRUE)
usethis::use_package("ggplot2", min_version = TRUE)

use_import_from("magrittr", "%>%")
use_import_from("rlang", ".data")


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


# NineSquares, the big one ------------------------------------------------

usethis::use_r("NineSquares")

document()
load_all()



## This goes into the examples of the function, but it's good to check that it works tbh
data <- data.frame(gene = paste0("Gene", 1:100),
                   LFC = rnorm(1000, mean = 5, sd = 2))
cutoffs <- NScutoff(data$LFC, scale = 1)



usethis::use_test("NineSquares")
test()

# checking behaviour of NSbasedf

data <- data.frame(
  gene = paste0("Gene", 1:10),
  num = c(2, 3, 4, 1, 5, 3, 2, 4, 5, 3),
  untreated_LFC = stats::rnorm(10, mean = 0, sd = 2),
  treated_LFC = stats::rnorm(10, mean = 0, sd = 2),
  untreated_pval = stats::runif(10, min = 0, max = 1),
  treated_pval = stats::runif(10, min = 0, max = 1)
)

NSbasedf(data, control = untreated_LFC,
         treament = treated_LFC,
         min_sgrna = 3)

check()
