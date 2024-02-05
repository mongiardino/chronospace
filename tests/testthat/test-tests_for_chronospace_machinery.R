test_that(desc = "tests for chronospace", code = {
  data("echinoid_dates")
  cspace<-chronospace(data_ages = echinoid_dates)

  expect_match(class(cspace), "chronospace")
  expect_equal(length(cspace), ncol(echinoid_dates$factors))

})
