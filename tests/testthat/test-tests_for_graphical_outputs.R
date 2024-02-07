test_that(desc = "tests for plot.chronospace", code = {
  data("echinoid_dates")
  cspace <- chronospace(data_ages = echinoid_dates)
  ordination <- plot.chronospace(cspace, output = "none")

  expect_match(class(ordination), "list")
  expect_equal(length(ordination), ncol(echinoid_dates$factors))
  for(i in 1:length(ordination)) expect_equal(names(ordination[[i]]), c("ordination", "PC_extremes"))
  for(i in 1:length(ordination)) expect_match(class(ordination[[i]]$ordination), "ggplot", all = FALSE)
  for(i in 1:length(ordination)) expect_equal(length(ordination[[i]]$PC_extremes), ncol(cspace[[i]]$ordination$x))
  for(i in 1:length(ordination)) {
    ord_i <- ordination[[i]]
    for(j in 1:length(ord_i$PC_extremes)) expect_match(class(ord_i$PC_extremes[[j]]), "ggtree", all = FALSE)
  }

})


test_that(desc = "tests for sensitive_nodes", code = {
  data("echinoid_dates")
  sensinodes <- sensitive_nodes(echinoid_dates, num_clades = 5, plot = FALSE)

  expect_equal(length(sensinodes), ncol(echinoid_dates$factors))
  for(i in 1:length(sensinodes)) expect_match(class(sensinodes[[i]]), "ggarrange", all = FALSE)

})


test_that(desc = "tests for ltt_sensitivity", code = {
  data("echinoid_dates")
  sensiltt <- ltt_sensitivity(echinoid_dates, plot = FALSE)

  expect_equal(length(sensiltt), ncol(echinoid_dates$factors))
  for(i in 1:length(sensiltt)) expect_match(class(sensiltt[[i]]), "ggplot", all = FALSE)

})


test_that(desc = "tests for specified_nodes", code = {
  data("echinoid_dates")
  specnode <- specified_node(echinoid_dates, tips = c('Brissus_obesus',
                                                      'Abatus_cordatus'), plot = TRUE)

  expect_equal(length(specnode), ncol(echinoid_dates$factors))
  for(i in 1:length(specnode)) expect_match(class(specnode[[i]]), "ggplot", all = FALSE)

})


