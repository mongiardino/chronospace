
#create chronospace-------------------------------------------------------------------

#' Create chronospace
#'
#' @description Compute the ordination maximizing variation in node ages using
#'   between-group PCA (one for each factor).
#'
#' @param data_ages A \code{"dataAges"} object created using [extract_ages()].
#' @param vartype Character, indicating the type of variation to be retained
#'   (\code{"total"} or \code{"non-redundant"}, see Details; not meaningful for
#'   factors with less than three levels).
#'
#' @details This function uses between-group PCA to find the set of axes
#'   maximizing variation in ages data between the groups of chronograms
#'   obtained through different methodological approaches. By default
#'   \code{vartype = "non-redundant"}, meaning bgPCA of each factor is performed
#'   using the variation left after removing the portion associated to all the
#'   other factors. If \code{vartype = "total"} (or if there is only one factor
#'   being assessed), bgPCA is performed over the raw variation in node ages.
#'
#' @return The total and non-redundant percentages of variation accounted for
#'   each factor are informed. An object of class \code{"chronospace"} is
#'   returned invisibly, with a list containing the ordination computed for each
#'   factor.
#'
#' @export
#'
#' @seealso \code{\link{plot.chronospace}}, \code{\link{sensitive_nodes}},
#'   \code{\link{ltt_sensitivity}}.
#'
#' @references
#'
#' @examples
#' #Load ages data
#' data("data_ages")
#'
#' #Create chronospace
#' cspace <- chronospace(data_ages)
#'
#' #Inspect object
#' cspace
chronospace <- function(data_ages, vartype = "non-redundant")  {

  #split data.frame 'data_ages' into ages and factors
  ages <- data_ages$ages
  factors <- data_ages$factors

  #factors names
  if(is.null(dim(factors))) factors <- data.frame(factors)
  facnames <- paste0("factor_", LETTERS[1:ncol(factors)])

  #create object for storing overall results, assign names
  results <- vector(mode = "list", length = ncol(factors))
  names(results) <- facnames

  #compute total variation
  totvar <- sum(apply(ages, 2, stats::var))

  #perform bgPCA for the crossing of all factors
  bgPCA0 <- bgprcomp(x = ages, groups = interaction(factors))
  fitted0 <- cbind(lm(bgPCA0$x ~ interaction(factors))$fitted)

  totexpvar <- sum(apply(fitted0, 2, stats::var))
  acc_tot <- 100 * (totexpvar / totvar)

  #perform bgPCA using each factor separately
  for(i in 1:ncol(factors)) {

    #perform bgPCA between groups defined by factor i over original variation
    bgPCA1 <- bgprcomp(x = ages, groups = factors[,i])
    fitted1 <- cbind(lm(bgPCA1$x ~ factors[,i])$fitted)

    #compute percentage of variation explained
    expvar <- sum(apply(fitted1, 2, stats::var))
    perc_tot <- 100 * (expvar / totvar)

    #report proportion of total variation explained
    if(ncol(factors) > 1) {
      catnames <- paste0('(', paste(unique(factors[,i]), collapse = '/'), ')')
      cat(paste('--- Results for', facnames[i], catnames, '---\n'))
    } else {
      cat('Results:\n')
    }

    cat(paste0('Proportion of total variation in node ages explained by ',
               facnames[i], ' = ',
               round(perc_tot, digits = 3),
               '%', '\n'))

    if(ncol(factors) > 1){
      #use bgPCA to compute an ordination that is residual to all factors but factor i
      bgPCA2.1 <- bgprcomp(x = ages, groups = interaction(factors[,-i]))
      resids2.1 <- bgPCA2.1$residuals

      #perform bgPCA between groups defined by factor i over residual variation
      bgPCA2.2 <- bgprcomp(x = resids2.1, groups = factors[,i])
      fitted2.2 <- cbind(lm(bgPCA2.2$x ~ factors[,i])$fitted)
      expvar2.2 <- sum(apply(fitted2.2, 2, stats::var))

      #compute percentage of non-redundant variation explained
      perc_nonred <- 100 * (expvar2.2 / totvar)

      #report proportion of non-redundant variation explained
      cat(paste0('Proportion of non-redundant variation in node ages explained by ',
                 facnames[i], ' = ',
                 round(perc_nonred, digits = 3),
                 '%', '\n\n'))
    } else {
      cat('(There is only one factor, non-redundant variation omitted)\n\n')
    }

    #select which bgPCA results are going to be used
    if(vartype == "total" | ncol(factors) == 1) {
      bgPCA <- bgPCA1
      perc <- perc_tot
    }

    if(vartype == "non-redundant" & ncol(factors) > 1) {
      bgPCA <- bgPCA2.2
      perc <- perc_nonred
    }


    #create table with percentages of variation explained by each axis
    vars <- apply(bgPCA$x, 2, stats::var) / totvar
    tab <- round(cbind(variance=vars, cummulative=cumsum(vars)), 5)

    #store bgPCA results, along with total variation and groups of factor i
    bgPCA$totvar <- totvar
    bgPCA$acc_tot <- acc_tot
    bgPCA$perc <- perc
    bgPCA$tab <- tab
    bgPCA$vartype <- vartype

    bgPCA$groups <- factors[,i]
    bgPCA$ages <- ages
    bgPCA$tree <- data_ages$topology
    results[[i]] <- bgPCA
  }

  class(results) <- "chronospace"

  return(invisible(results))
}

# print chronospace object -------------------------------------------------------

#' Print \code{"dataAges"} objects
#'
#' @param x A \code{"chronospace"} object
#'
#' @return Information on percentages of age variation explained by the included
#'   factors
#'
#' @export
#'
#' @examples
#' #Load ages data
#' data("data_ages")
#'
#' #Create chronospace
#' cspace <- chronospace(data_ages)
#'
#' #Inspect object
#' print(cspace)
print.chronospace <- function(x) {
  cat(paste0("Percentage of variation accounted by all factors (including interactions): ",
             round(x[[1]]$acc_tot, 3), "%", "\n"))
  for(i in 1:length(x)) {
    cat(paste0("Percentage of ", x[[i]]$vartype, " variation accounted by ", names(x)[i], " : ",
               round(x[[i]]$perc, 3), "%", "\n"))
    tab <- as.data.frame(x[[i]]$tab)
    rownames(tab) <- paste0("bgPC",1:nrow(tab))
    print(tab)
  }

}



# internal between-group PCA function ---------------------------------------------
bgprcomp <- function(x, groups) {

  grandmean <- colMeans(x)
  x_centered <- scale(x, scale = F, center = T)
  x_gmeans <- apply(X = x_centered, MARGIN = 2, FUN = tapply, groups, mean)

  V_g <- stats::cov(x_gmeans)
  eig <- eigen(V_g)

  scores <- x_centered %*% eig$vectors
  scores <- cbind(scores[,1:(nlevels(groups) - 1)])
  rotation <- eig$vectors

  preds <- scores %*% t(rotation[,1:ncol(scores)])
  resids <- x - preds

  return(list(x = scores, residuals = resids, rotation = rotation,
              values = eig$values, center = grandmean, gmeans = x_gmeans))
}


# internal reverse PCA function -----------------------------------------------------
revPCA <- function(scores, vectors, center) { t(t(scores %*% t(vectors)) + center) }


#translate ages into trees-------------------------------------------------------------------

reconstruct_blen <- function(clades, tree, plus, mean, minus) {

  #check number of descendants stemming from each node
  clade_size <- unlist(lapply(clades, length))

  #setup trees that will have mean branch lengths, mean+sdev and mean-sdev
  #trees
  tree_mean <- tree_plus <- tree_minus <- tree

  #loop through clades from smallest to biggest (i.e., up the tree)
  for(k in 2:max(clade_size)) {
    #which nodes have the number of descendants
    which_clades <- which(clade_size == k)
    if(length(which_clades) > 0) {
      for(l in 1:length(which_clades)) {
        #which node are we talking about
        node_to_change <- ape::getMRCA(tree, unlist(clades[which_clades[l]]))

        #get node ages for this node
        dif_minus <- minus[,which_clades[l]]
        dif_mean <- mean[which_clades[l],]
        dif_plus <- plus[,which_clades[l]]

        #if the clade is a cherry (i.e., 2 descendants)
        if(k == 2) {
          #get branches descending to both tips and assign them their
          #correct branches (which is == to the node age)
          branches_to_descendants <- which(tree$edge[,1] == node_to_change)
          tree_minus$edge.length[branches_to_descendants] <- dif_minus
          tree_mean$edge.length[branches_to_descendants] <- dif_mean
          tree_plus$edge.length[branches_to_descendants] <- dif_plus

          #if it is not a cherry
        } else {
          #get nodes of direct descendant
          nodes_of_descendants <- tree$edge[,2][which(tree$edge[,1] ==
                                                        node_to_change)]

          #if any descendant is a tip do as above, assign branch length ==
          #node age
          if(any(nodes_of_descendants %in% 1:length(tree$tip.label))) {
            singletons <- nodes_of_descendants[which(nodes_of_descendants %in%
                                                       1:length(tree$tip.label))]
            branches_to_singletons <- which(tree$edge[,2] == singletons)
            tree_minus$edge.length[branches_to_singletons] <- dif_minus
            tree_mean$edge.length[branches_to_singletons] <- dif_mean
            tree_plus$edge.length[branches_to_singletons] <- dif_plus

            #remove it from descendants as its branch length is already
            #set
            nodes_of_descendants <- nodes_of_descendants[-which(nodes_of_descendants ==
                                                                  singletons)]
          }

          #for descendant clades do the following
          for(m in 1:length(nodes_of_descendants)) {
            #obtain all descendants
            tips <- unlist(phangorn::Descendants(tree, nodes_of_descendants[m],
                                                 type = 'tips'))

            #remove from the age the age of the descendant node, which is
            #already set up correctly as the loop goes from smaller to
            #larger clades
            tree_minus$edge.length[which(tree_minus$edge[,2] == nodes_of_descendants[m])] <-
              dif_minus - ape::dist.nodes(tree_minus)[tips[1], nodes_of_descendants[m]]
            tree_mean$edge.length[which(tree_mean$edge[,2] == nodes_of_descendants[m])] <-
              dif_mean - ape::dist.nodes(tree_mean)[tips[1], nodes_of_descendants[m]]
            tree_plus$edge.length[which(tree_plus$edge[,2] == nodes_of_descendants[m])] <-
              dif_plus - ape::dist.nodes(tree_plus)[tips[1], nodes_of_descendants[m]]
          }
        }
      }
    }
  }
  return(list(tree_plus = tree_plus, tree_mean = tree_mean, tree_minus = tree_minus))
}
