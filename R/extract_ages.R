#load input data ----------------------------------------------------------

#' Extract and organize ages data
#'
#' @description Extract and organize ages from files containing the posterior
#'   distributions of time-calibrated phylogenetic trees (i.e. chronograms)
#'   obtained through different methodological decisions.
#'
#' @importFrom magrittr %>%
#'
#' @param path A character string specifying the path to the folder containing
#'   (only) the input files.
#' @param type  A list of \code{f} vectors (one for each factor being tested)
#'   specifying the group to which the chronograms from each file will be
#'   assigned to.
#' @param sample The fixed number of trees to retain from each file.
#'
#' @details This function imports and prepare the ages data such that they are
#'   suitable for analysis using chronospace. The input files are a series of
#'   files in Newick format containing time-calibrated trees taken from the
#'   posterior distributions of separate Bayesian inferences of the same dataset
#'   (can be either different runs of the same analysis or separate analyses).
#'   These should differ in specific aspects of their initial setting whose
#'   effects are of interest (e.g. the model of evolution used and/or the method
#'   employed to subsample genes) and be stored in an exclusive folder. The
#'   topology is assumed constrained, such that chronograms differ only in the
#'   branch lengths (and therefore, in their inferred node ages).
#'
#' @return A \code{t x (n + f)} matrix in which each row is one of \code{t}
#'   trees, and each of the first \code{n} clumns represents a node in the fixed
#'   topology. The last \code{f} columns indicate classification of each tree to
#'   each factor.
#'
#' @export
#'
#' @seealso \code{\link{chronospace}}
#'
#' @examples
extract_ages <- function(path = NA, type, sample) {

  #Obtain the names of all files in the specified directory (or the working
  #directory otherwise). These should all be newick tree files corresponding to
  #Bayesian posterior distributions of time-calibrated analyses with constrained
  #topology.
  if(is.na(path)) {
    files <- list.files()
  } else {
    files <- list.files(path = path)
  }

  #check that tree files and factors provided in 'types' match correctly
  cat("Check that labels are assigned correctly to the input files.\n",
      "If there is an error, modify the order of factors in 'type'\n",
      "or the name of input files, for the two to match.\n\n", sep = '')

  for(i in 1:length(files)) {
    to_print <- paste0('file = ', files[i], ' | type = ',
                       paste0(unname(as.vector(as.data.frame(type)[i,])),
                              collapse = ' - '))
    cat(to_print, '\n')
  }

  #loop through the tree files, load them and subsample them to the number
  #specified in 'sample'
  for(i in 1:length(files)) {
    if(is.na(path)) {
      trees <- ape::read.tree(paste0(getwd(), '/', files[i]))
    } else {
      trees <- ape::read.tree(paste0(path, '/', files[i]))
    }
    if(!is.na(sample)) {
      trees <- trees[sample(1:length(trees), sample)]
    }

    if(i == 1) {
      all_trees <- trees
    } else {
      all_trees <- c(all_trees, trees)
    }
  }

  #extract the number of nodes in the tree and the taxonomic composition of each
  #node
  tree <- all_trees[[1]]
  clades <- list()
  for(i in 1:tree$Nnode) {
    clades[i] <- list(tree$tip.label[unlist(phangorn::Descendants(tree, length(tree$tip.label)+i, type = 'tips'))])
  }

  #build the matrix that will contain node ages (columns) for each tree (rows)
  ages <- matrix(0, ncol = length(clades), nrow = length(all_trees))

  #assign values from the oldest node (i.e., root)
  root <- sapply(all_trees, function(x) max(phytools::nodeHeights(x)))
  ages[,1] <- root

  #assign values to all other nodes
  for(i in 2:ncol(ages)) {
    #check node number associated with clade i
    node <- sapply(all_trees, phytools::findMRCA, clades[[i]])

    #if the clade is always associated with the same number this runs fast
    if(length(unique(node)) == 1) {
      node <- node[1]
      pos <- which(all_trees[[1]]$edge[,2] == node)
      ages[,i] <- root - sapply(all_trees, function(x) phytools::nodeHeights(x)[pos,2])

      #otherwise this takes a while, but the correct matching is confirmed such
      #that dates from the same clade are placed in the same column regardless
      #of its node number
    } else {
      for(j in 1:length(all_trees)) {
        ages[j,i] <- root[j] -
          phytools::nodeHeights(all_trees[[j]])[which(all_trees[[j]]$edge[,2] == node[j]),2]
      }
    }
  }

  #change to data.frame, set column names and add types of runs as factors
  data_ages <- data.frame(ages)
  colnames(data_ages) <- paste0('clade_', 1:ncol(data_ages))

  types_of_runs <- data.frame(matrix(NA, nrow = nrow(ages), ncol = length(type)))
  for(i in 1:length(type)) {
    types_of_runs[,i] <- rep(type[[i]], each = nrow(ages)/length(type[[i]]))
  }

  colnames(types_of_runs) <- paste('factor', LETTERS[1:ncol(types_of_runs)],
                                   sep = '_')

  data_ages <- cbind(data_ages, types_of_runs)
  data_ages <- data_ages %>% dplyr::mutate_if(sapply(data_ages, is.character), as.factor)

  #export
  return(data_ages)
}

#load input data ----------------------------------------------------------

#' Extract and organize ages data
#'
#' @description Extract and organize ages from files containing the posterior
#'   distributions of time-calibrated phylogenetic trees (i.e. chronograms)
#'   obtained through different methodological decisions.
#'
#' @importFrom magrittr %>%
#'
#' @param path A character string specifying the path to the folder containing
#'   (only) the input files.
#' @param type  A list of \code{f} vectors (one for each factor being tested)
#'   specifying the group to which the chronograms from each file will be
#'   assigned to.
#' @param sample The fixed number of trees to retain from each file.
#'
#' @details This function imports and prepare the ages data such that they are
#'   suitable for analysis using chronospace. The input files are a series of
#'   files in Newick format containing time-calibrated trees taken from the
#'   posterior distributions of separate Bayesian inferences of the same dataset
#'   (can be either different runs of the same analysis or separate analyses).
#'   These should differ in specific aspects of their initial setting whose
#'   effects are of interest (e.g. the model of evolution used and/or the method
#'   employed to subsample genes) and be stored in an exclusive folder. The
#'   topology is assumed constrained, such that chronograms differ only in the
#'   branch lengths (and therefore, in their inferred node ages).
#'
#' @return An object of class \code{"dataAges"} with a list containing 1) the
#'   \code{t x n} ages matrix (where each row is one of \code{t} trees and each
#'   column one of \code{n} nodes), 2) the \code{f} associated factors, and 3)
#'   the fixed phylogenetic topology.
#'
#' @export
#'
#' @seealso \code{\link{chronospace}}
#'
#' @examples
extract_ages2 <- function(path = NA, type, sample) {

  #Obtain the names of all files in the specified directory (or the working
  #directory otherwise). These should all be newick tree files corresponding to
  #Bayesian posterior distributions of time-calibrated analyses with constrained
  #topology.
  if(is.na(path)) {
    files <- list.files()
  } else {
    files <- list.files(path = path)
  }

  #check that tree files and factors provided in 'types' match correctly
  cat("Check that labels are assigned correctly to the input files.\n",
      "If there is an error, modify the order of factors in 'type'\n",
      "or the name of input files, for the two to match.\n\n", sep = '')

  for(i in 1:length(files)) {
    to_print <- paste0('file = ', files[i], ' | type = ',
                       paste0(unname(as.vector(as.data.frame(type)[i,])),
                              collapse = ' - '))
    cat(to_print, '\n')
  }

  #loop through the tree files, load them and subsample them to the number
  #specified in 'sample'
  for(i in 1:length(files)) {
    if(is.na(path)) {
      trees <- ape::read.tree(paste0(getwd(), '/', files[i]))
    } else {
      trees <- ape::read.tree(paste0(path, '/', files[i]))
    }
    if(!is.na(sample)) {
      trees <- trees[sample(1:length(trees), sample)]
    }

    if(i == 1) {
      all_trees <- trees
    } else {
      all_trees <- c(all_trees, trees)
    }
  }

  #extract the number of nodes in the tree and the taxonomic composition of each
  #node
  tree <- all_trees[[1]]
  clades <- list()
  for(i in 1:tree$Nnode) {
    clades[i] <- list(tree$tip.label[unlist(phangorn::Descendants(tree, length(tree$tip.label)+i, type = 'tips'))])
  }

  #build the matrix that will contain node ages (columns) for each tree (rows)
  ages <- matrix(0, ncol = length(clades), nrow = length(all_trees))

  #assign values from the oldest node (i.e., root)
  root <- sapply(all_trees, function(x) max(phytools::nodeHeights(x)))
  ages[,1] <- root

  #assign values to all other nodes
  for(i in 2:ncol(ages)) {
    #check node number associated with clade i
    node <- sapply(all_trees, phytools::findMRCA, clades[[i]])

    #if the clade is always associated with the same number this runs fast
    if(length(unique(node)) == 1) {
      node <- node[1]
      pos <- which(all_trees[[1]]$edge[,2] == node)
      ages[,i] <- root - sapply(all_trees, function(x) phytools::nodeHeights(x)[pos,2])

      #otherwise this takes a while, but the correct matching is confirmed such
      #that dates from the same clade are placed in the same column regardless
      #of its node number
    } else {
      for(j in 1:length(all_trees)) {
        ages[j,i] <- root[j] -
          phytools::nodeHeights(all_trees[[j]])[which(all_trees[[j]]$edge[,2] == node[j]),2]
      }
    }
  }

  #change to data.frame, set column names and add types of runs as factors
  data_ages <- data.frame(ages)
  colnames(data_ages) <- paste0('clade_', 1:ncol(data_ages))

  #set information about factors
  types_of_runs <- data.frame(matrix(NA, nrow = nrow(ages), ncol = length(type)))
  for(i in 1:length(type)) {
    types_of_runs[,i] <- rep(type[[i]], each = nrow(ages)/length(type[[i]]))
  }

  colnames(types_of_runs) <- paste('factor', LETTERS[1:ncol(types_of_runs)],
                                   sep = '_')

  types_of_runs <- types_of_runs %>% dplyr::mutate_if(sapply(types_of_runs, is.character), as.factor)

  #strip tree from branch lengths
  tree$edge.length<-NULL

  #put results together
  results<-list(ages=data_ages, factors=types_of_runs, topology=tree)

  #export
  class(results)<-c("list", "dataAges")
  return(results)
}

# print ages data --------------------------------------------------

#' Print dataAges objects
#'
#' @param x a dataAges object.
#'
#' @return
#'
#' @export
#'
#' @examples
print.dataAges<-function(x) {

  ages<-x$ages
  facs<-x$factors

  infoages<-paste0("Data from ", nrow(ages),
                   " trees with ", ncol(ages),
                   " internal nodes (see $ages and $topology),")

  if(ncol(facs)==1) {
    infofacs<-paste0("Obtained using 1 methodological pathway (see $factors).")
  } else {
    infofacs<-paste0("Obtained using ", ncol(facs), " methodological pathways (see $factors).")
  }

  cat(infoages, infofacs, sep="\n")

}

