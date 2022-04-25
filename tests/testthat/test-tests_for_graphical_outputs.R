test_that(desc = "tests for plot.chronospace", code = {
  data("data_ages")
  cspace <- chronospace(data_ages = data_ages)
  ordination <- plot(obj=cspace)

  expect_match(class(ordination), "list")
  expect_equal(length(ordination), length(cspace))
  for(i in 1:length(ordination)) expect_equal(names(ordination[[i]]), c("ordination", "PC_extremes"))
  for(i in 1:length(ordination)) expect_match(class(ordination[[i]]$ordination), "ggplot", all = FALSE)
  for(i in 1:length(ordination)) expect_equal(length(ordination[[i]]$PC_extremes), ncol(cspace[[i]]$x))
  for(i in 1:length(ordination)) {
    ord_i <- ordination[[i]]
    for(j in 1:length(ord_i$PC_extremes)) expect_match(class(ord_i$PC_extremes[[j]]), "ggtree", all = FALSE)
  }

})


test_that(desc = "tests for sensitive_nodes", code = {
  data("data_ages")
  cspace <- chronospace(data_ages = data_ages)
  sensinodes <-sensitive_nodes(obj = cspace, chosen_clades = 5)

  expect_equal(length(sensinodes), length(cspace))
  for(i in 1:length(sensinodes)) expect_match(class(sensinodes[[i]]), "ggarrange", all = FALSE)

})


test_that(desc = "tests for ltt_sensitivity", code = {
  data("data_ages")
  sensiltt<-ltt_sensitivity(data_ages = data_ages, average = "mean")

  expect_equal(length(sensiltt), ncol(data_ages$factors))
  for(i in 1:length(sensiltt)) expect_match(class(sensiltt[[i]]), "ggplot", all = FALSE)

})
