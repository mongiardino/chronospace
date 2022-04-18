#load input data ----------------------------------------------------------

#' Extract and organize ages data
#'
#' @description Extract and organize ages from files containing the posterior distributions of time-calibrated
#' phylogenetic trees (i.e. chronograms) obtained through different methodological decisions.
#'
#' @importFrom magrittr %>%
#'
#' @param path A character string specifying the path to the folder containing (only) the input files.
#' @param type  A list of g vectors (one for each factor being tested) specifying the group to which the chronograms
#'  from each file will be assigned to.
#' @param sample The fixed number of trees to retain from each file.
#'
#' @details This function imports and prepare the ages data such that they are suitable for analysis using chronospace.
#' The input files are a series of files in Newick format containing time-calibrated trees taken from the posterior
#' distributions of separate Bayesian inferences of the same dataset (can be either different runs of the same analysis
#' or separate analyses). These should differ in specific aspects of their initial setting whose effects are of interest
#' (e.g. the model of evolution used and/or the method employed to subsample genes) and be stored in an exclusive folder.
#' The topology is assumed constrained, such that chronograms differ only in the branch lengths (and therefore, in their
#' inferred node ages).
#'
#' @return A t x (n + g) matrix in which each row is one of t trees, and each of the first n clumns reprresents a node
#' in the fixed topology. The last g columns indicate classification of each tree to each of the methodological groups.
#' @export
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
      "or the name of input files for the two to match.\n\n", sep = '')

  for(i in 1:length(files)) {
    to_print <- paste0('file = ', files[i], ' | type = ',
                       paste0(unname(as.vector(as.data.frame(type)[i,])),
                              collapse = ' - '))
    cat(to_print, '\n')
  }

  #loop through the tree files, load them and subsample them to the number
  #specified in 'sample'
  for(i in 1:length(files)) {
    trees <- ape::read.tree(paste0(getwd(), '/', files[i]))
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
    node <- sapply(all_trees, findMRCA, clades[[i]])

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
