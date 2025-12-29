
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


# Startup -----------------------------------------------------------------

# Should come preloaded from the .Rprofile thing

library(devtools)


# First function: %nin% ---------------------------------------------------

usethis::use_r("nin")

# here I added %nin% code to the R file and saved it.

load_all() # should make the %nin% function available

# we can check it works after loading
c("A", "B", "C", "D", "E", "F") %nin% c("B", "A", "D", "E")

exists("%nin%", where = globalenv(), inherits = FALSE)
# has to throw FALSE because we're in the development environment, not global


check()
