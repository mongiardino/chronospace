test_that(desc = "tests for chronospace", code = {
  data("data_ages")
  cspace<-chronospace(data_ages = data_ages)

  expect_match(class(cspace), "chronospace")
  expect_equal(length(cspace), ncol(data_ages$factors))

})
