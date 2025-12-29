test_that("%nin% throws intended boolean", {
  expect_equal(
    c("A", "B", "C", "D", "E", "F") %nin% c("B", "A", "D", "E"),
    c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE)
  )
  expect_all_false("samba" %nin% c("salsa", "samba", "son"))
  expect_all_false(c(1, 3, 5, 7) %nin% 1:10)
})

test_that("%nin% is case sensitive", {
  expect_all_true(c("a", "b", "c", "d", "e", "f") %nin% c("B", "A", "D", "E"))
  expect_all_true("SaLsA" %nin% c("salsa", "samba", "son"))
  expect_all_true("salsa" %nin% c("Salsa", "sALSa", "saLsa"))
})


test_that("%nin% handles empty vectors", {
  expect_equal(character(0) %nin% c("a", "b"), logical(0))
})

test_that("%nin% handles NA correctly", {
  expect_equal(c("a", NA) %nin% c("a", "b"), c(FALSE, TRUE))
})
