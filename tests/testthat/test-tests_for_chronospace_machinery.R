test_that(desc = "tests for chronospace", code = {
  type <- list(c('clock', 'clock', 'random', 'random', 'signal', 'signal'),
               c('CATGTR', 'GTR', 'CATGTR', 'GTR', 'CATGTR', 'GTR'))
  data <- extract_ages2(path="C:/Users/pablo/Desktop/Trabajo/Investigacion/Proyectos/chronospace_local/data",
                        type = type, sample = 500)
  cspace<-chronospace2(data_ages = data, vartype = "non-redundant")

  expect_match(class(cspace), "chronospace")
  expect_equal(length(cspace), ncol(data$factors))

  rm(list=c("data","cspace"))
})
