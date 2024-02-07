
#summarize chronospace----------------------------------------------------------

#' Summarize chronospace
#'
#' @description Compute percentages of variation in node ages accounted by each
#'   factor, and extract axes of variation maximizing separation between factor
#'   levels using between-group PCA.
#'
#' @param data_ages A \code{"nodeAges"} object created using [extract_ages()].
#'
#' @details This function summarizes variation in node ages in two ways. First,
#'   a between-groups PCA is performed for each individual factor, extracting
#'   ordination axes that maximize separation between the levels of that factor,
#'   which can be visualized using \code{\link{plot.chronospace}}. Second, a Sum
#'   of squares approaches is used to quantify variation accounted by each
#'   factor, both from the raw total variation and from the bgPC axes computed
#'   for each individual factor.
#'
#' @return Percentages of variation in node ages accounted by the factors
#'   considered, both in total and individually, are printed. An object of
#'   class \code{"chronospace"} is returned invisibly, with a list containing
#'   the ordination computed for each factor.
#'
#' @export
#'
#' @seealso \code{\link{plot.chronospace}}
#'
#' @references
#'
#' @examples
#' #Load ages data
#' data("echinoid_dates")
#'
#' #Create chronospace
#' cspace <- chronospace(echinoid_dates)
#'
#' #print table
#' cspace
#'
#' #access ordination and SSQ for factor A
#' cspace$factor_A$ordination
chronospace <- function(data_ages) {

  #split data.frame 'data_ages' into ages and factors
  ages <- data_ages$ages
  factors <- data_ages$factors

  #extract factors names
  facnames <- colnames(factors)

  #create object for storing overall results, assign names
  results <- vector(mode = "list", length = ncol(factors))
  names(results) <- facnames

  #compute total variation (sum trace of covariance matrix)
  totvar <- sum(apply(ages, 2, stats::var))


  #decompose total variation in node ages. For simplification, interactions are
  #assumed to be absent.
  form <- formula(paste0("as.matrix(ages)[,i] ~", paste(facnames, collapse = " + ")))
  n <- nrow(ages)

  SS <- NULL
  for(i in 1:ncol(ages)) {
    model <- lm(form, data = data.frame(factors, ages))
    ss <- anova(model)$`Sum Sq`/ (n - 1)
    SS <- rbind(SS, ss)
  }
  colnames(SS) <- paste0(row.names(anova(model)), sep = " (%)")
  colnames(SS)[colnames(SS) == "Residuals (%)"] <- "Unaccounted (%)"
  SS_perc <- rbind(round(100 * colSums(SS)/totvar, 5))
  rownames(SS_perc) <- "Total variation"

  print(SS_perc)
  cat("_________________________________________________________________________________\n")

  vartable <- data.frame("All", SS_perc)
  colnames(vartable) <- c("Grouping_factor", colnames(SS))

  #perform bgPCA using each factor separately
  for(i in 1:ncol(factors)) {

    #perform bgPCA between groups defined by factor i over original variation
    bgPCA <- bgprcomp(x = ages, groups = factors[,i])

    subtotvar <- sum(apply(bgPCA$x, 2, var))
    form <- formula(paste0("bgPCA$x[,j] ~", paste(facnames, collapse = " + ")))

    subSS <- NULL
    for(j in 1:ncol(bgPCA$x)) {
      model <- stats::lm(form, data = data.frame(factors))
      ss <- stats::anova(model)$`Sum Sq`/ (n - 1)
      subSS <- rbind(subSS, ss)
    }

    colnames(subSS) <- paste0(row.names(stats::anova(model)), sep = " (%)")
    colnames(subSS)[colnames(subSS) == "Residuals (%)"] <- "Unaccounted (%)"
    subSS_perc <- round(100 * subSS / totvar, 5)
    rownames(subSS_perc) <- paste0("bgPC", 1:ncol(bgPCA$x), "(",
                                   round(
                                     100 * apply(bgPCA$x, 2, var) /
                                       totvar, 2), "%)")

    #report proportion of total variation explained
    if(ncol(factors) > 1) {
      catnames <- paste0('(', paste(unique(factors[,i]), collapse = '/'), ')')
      cat(paste('--- Results for', facnames[i], catnames, '---\n'))
    } else {
      cat('Results:\n')
    }

    #print
    print(subSS_perc)
    cat("---------------------------------------------------------------------------------\n")

    subvartable <- data.frame(facnames[i], subSS_perc)
    colnames(subvartable) <- c("Grouping_factor", colnames(subSS))


    #store bgPCA results, along with total variation and groups of factor i

    result <- list()
    result$ordination <- bgPCA
    result$ssq <- list(totvar = totvar, vartable = subvartable)
    result$data <- list(ages = ages, groups = factors[,i],
                        tree = data_ages$topology)

    results[[i]] <- result

  }

  results$Total_vartable <- vartable

  cat(" * All percentages are relative to the total amount of variation in node ages\n")

  class(results) <- "chronospace"

  return(invisible(results))
}


# print chronospace object -------------------------------------------------------

#' Print \code{"chronospace"} objects
#'
#' @param x A \code{"chronospace"} object.
#'
#' @return Reports percentages in node ages variation explained by the factors
#'   included.
#'
#' @export
#'
#' @examples
#' #Load ages data
#' data("echinoid_dates")
#'
#' #Create chronospace
#' cspace <- chronospace(echinoid_dates)
#'
#' #Inspect object
#' print(cspace)
print.chronospace <- function(x) {
  cat("_____________________________________________________________________________________\n")
  print(x$Total_vartable)
  cat("_____________________________________________________________________________________\n")
  x <- x[names(x) != "Total_vartable"]
  for(i in seq_len(length(x))) {
    print(x[[i]]$ssq$vartable)
    cat("-------------------------------------------------------------------------------------\n")
  }
  cat(" * All percentages are relative to the total amount of variation in node ages\n")
}


# internal between-group PCA function ---------------------------------------------
bgprcomp <- function(x, groups) {

  grandmean <- colMeans(x)
  x_centered <- scale(x, scale = FALSE, center = TRUE)
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
