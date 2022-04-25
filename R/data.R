#' Sample of chronograms from the Echinodermata phylogeny
#'
#' @description Ages data of 3000 chronograms extracted from six posterior
#'   distributions of time-calibrated trees, created using [extract_ages()].
#'
#' @format a \code{"dataAges"} object with a list containing
#' \describe{
#'   \item{$ages}{a \code{t x n} data.frame with ages data (where each row is
#'   one of \code{t} trees and each column one of \code{n} nodes.}
#'   \item{$factors}{a \code{t x f} data.frame specifying the classification
#'   of each of the \code{t} trees into \code{f} factors.}
#'   \item{$topology}{the fixed phylogenetic topology.}
#'   }
#'
#' @details A sample of chronograms extracted from six posterior distributions
#'   of time-calibrated trees (500 trees each), obtained using PhyloBayes
#'   (Lartillot et al. 2013). These were all run using the same contrained
#'   topology amd three different sets of 100 genes subsampled from a larger
#'   phylogenomic dataset based on their level of clock-likeness or phylogenetic
#'   signal, or otherwise at random. Each subsample was also run under two
#'   models of molecular evolution, the site-homogeneous GTR+G and the
#'   site-heterogeneous CAT+GTR+G (Lartillot & Philippe 2004).
#'
#' @references
"data_ages"