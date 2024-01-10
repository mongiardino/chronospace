#load input data ----------------------------------------------------------

#' Import and format node ages data
#'
#' @description Extract and organize ages from files containing the posterior
#'   distributions of time-calibrated phylogenetic trees (chronograms) obtained
#'   through different methodological decisions.
#'
#' @importFrom magrittr %>%
#'
#' @param path An optional character string specifying the path to the folder
#'   containing (only) the input files. If \code{NULL}, the current working
#'   directory is used instead.
#' @param type A list of \code{f} vectors (one for each factor being tested)
#'   specifying the group to which the chronograms from each file will be
#'   assigned to.
#' @param sample Numeric; the fixed number of trees to retain from each file.
#' @param facnames An optional character vector indicating the names of the
#'   factors.
#'
#' @details This function imports and prepare the ages data such that they are
#'   suitable for analysis using the chronospace suite of functions. The input
#'   files are a series of files in Newick format containing time-calibrated
#'   trees taken from the posterior distributions of separate Bayesian
#'   inferences of the same data set (can be either different runs of the same
#'   analysis or separate analyses). These should differ in specific aspects of
#'   their initial setting whose effects are of interest (e.g. the model of
#'   evolution used and/or the method employed to subsample genes), and be
#'   stored in an exclusive folder. The topology is assumed to be constrained
#'   such that chronograms differ only in branch length (and therefore in their
#'   inferred node ages).
#'
#' @return An object of class \code{"nodeAges"} with a list containing 1) the
#'   \code{t x n} ages matrix (where each row is one of \code{t} trees and each
#'   column one of \code{n} nodes), 2) the \code{f} associated factors, and 3)
#'   the fixed phylogenetic topology.
#'
#' @export
#'
#' @seealso \code{\link{chronospace}}, \code{\link{sensitive_nodes}},
#'   \code{\link{ltt_sensitivity}}
#'
#' @examples
#' #Create temporal directory for the files
#' temp0 <- tempdir()
#' dir.create(paste(temp0, "files", sep = "/"))
#' temp <- paste(temp0, "files", sep = "/")
#'
#' #Declare and download files from the internet (this might take a minute)
#' url <- "https://raw.githubusercontent.com/mongiardino/chronospaces_eLife/main/example_files/"
#' files <- c("clockCATGTR_ln_sample.datedist",
#'            "clockGTR_ln_sample.datedist",
#'            "randomCATGTR_ln_sample.datedist",
#'            "randomGTR_ln_sample.datedist",
#'            "signalCATGTR_ln_sample.datedist",
#'            "signalGTR_ln_sample.datedist")
#' for(i in 1:length(files)) download.file(paste0(url, files[i]), paste(temp, files[i], sep = "/"))
#'
#' #Check files names, compare against type order below
#' list.files(temp)
#'
#' #Set type of runs and number of chronograms to be retained
#' type <- list(c('clock', 'clock', 'random', 'random', 'signal', 'signal'),
#'              c('CATGTR', 'GTR', 'CATGTR', 'GTR', 'CATGTR', 'GTR'))
#' sample <- 500
#'
#' #Import data to R (this might take a minute)
#' data <- extract_ages(path = temp, type = type, sample = sample)
#' data
extract_ages <- function(path = NULL, type, sample, facnames = NULL) {


  #Obtain the names of all files in the specified directory (or the working
  #directory otherwise). These should all be Newick tree files corresponding to
  #Bayesian posterior distributions of time-calibrated analyses with constrained
  #topology.
  if(is.null(path)) {
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

    if(is.null(path)) {
      trees <- ape::read.tree(paste0(getwd(), '/', files[i]))
    } else {
      trees <- ape::read.tree(paste0(path, '/', files[i]))
    }
    if(!is.null(sample)) {
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
    clades[i] <- list(tree$tip.label[unlist(phangorn::Descendants(tree, length(tree$tip.label) + i, type = 'tips'))])
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
    types_of_runs[,i] <- rep(type[[i]], each = nrow(ages) / length(type[[i]]))
  }

  if(is.null(facnames)) {
    colnames(types_of_runs) <- paste('factor', LETTERS[1:ncol(types_of_runs)],
                                     sep = '_')
  } else {
    colnames(types_of_runs) <- facnames
  }

  types_of_runs <- types_of_runs %>% dplyr::mutate_if(sapply(types_of_runs, is.character),
                                                      as.factor)

  #strip tree from branch lengths
  tree$edge.length <- NULL

  #put results together
  results<-list(ages = data_ages, factors = types_of_runs, topology = tree)

  #export
  class(results) <- c("list", "nodeAges")
  return(results)
}

# print ages data --------------------------------------------------

#' Print \code{"nodeAges"} objects
#'
#' @param x a \code{"nodeAges"} object.
#'
#' @return Reports the number of trees, factors, and general structure of
#'   the object.
#'
#' @export
#'
#' @examples
#' data("data_ages")
#' print(data_ages)
print.nodeAges <- function(x) {

  ages <- x$ages
  npathw <- prod(apply(x$factors, 2, \(x) {length(unique(x))}))

  infoages <- paste0("Data from ", nrow(ages),
                     " trees with ", ncol(ages),
                     " internal nodes (see $ages and $topology),")

  infofacs <- paste0("Obtained using ", npathw, " methodological pathways (see $factors).")
  cat(infoages, infofacs, sep = "\n")

}

