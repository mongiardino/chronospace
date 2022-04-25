test_that(desc = "tests for extract_ages", code = {
  type <- list(c('clock', 'clock', 'random', 'random', 'signal', 'signal'),
               c('CATGTR', 'GTR', 'CATGTR', 'GTR', 'CATGTR', 'GTR'))
  data <- extract_ages2(path="C:/Users/pablo/Desktop/Trabajo/Investigacion/Proyectos/chronospace_local/data",
                        type = type, sample = 500)

  expect_match(class(data), "dataAges", all = FALSE)
  expect_match(class(data), "list", all = FALSE)
  expect_equal(data$topology$Nnode, ncol(data$ages))

  rm(data)
})
