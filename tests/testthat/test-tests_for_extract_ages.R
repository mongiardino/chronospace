test_that(desc = "tests for extract_ages", code = {
  temp0 <- tempdir()
  dir.create(paste(temp0, "files", sep="\\"))
  temp <- paste(temp0, "files", sep="\\")

  url <- "https://raw.githubusercontent.com/mongiardino/chronospace/main/example_files/"
  files <- c("clockCATGTR_ln_sample.datedist",
             "clockGTR_ln_sample.datedist",
             "randomCATGTR_ln_sample.datedist",
             "randomGTR_ln_sample.datedist",
             "signalCATGTR_ln_sample.datedist",
             "signalGTR_ln_sample.datedist")
  for(i in 1:length(files)) utils::download.file(paste0(url, files[i]), paste(temp, files[i], sep="\\"))

  type <- list(c('clock', 'clock', 'random', 'random', 'signal', 'signal'),
               c('CATGTR', 'GTR', 'CATGTR', 'GTR', 'CATGTR', 'GTR'))
  sample <- 500

  data <- extract_ages(path = temp, type = type, sample = sample)

  expect_match(class(data), "dataAges", all = FALSE)
  expect_match(class(data), "list", all = FALSE)
  expect_equal(data$topology$Nnode, ncol(data$ages))

})




