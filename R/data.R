#' Sample of chronograms from the Echinodermata phylogeny
#'
#' @description Ages data of 3000 chronograms extracted from six posterior
#'   distributions of time-calibrated trees, created using [extract_ages()].
#'
#' @format a \code{"nodeAges"} object with a list containing
#' \describe{
#'   \item{$ages}{: a \code{t x n} data frame containing node ages data (where
#'   each row is one of \code{t} trees and each column one of \code{n} nodes.}
#'   \item{$factors}{: a \code{t x f} data frame containing the classification
#'   of each of the \code{t} trees into \code{f} factors.}
#'   \item{$topology}{: the fixed phylogenetic topology.}
#'   }
#'
#' @details A sample of chronograms extracted from six posterior distributions
#'   of time-calibrated trees (500 trees each), obtained using PhyloBayes
#'   (Lartillot et al. 2013). These were all run using the same constrained
#'   topology and three different sets of 100 genes subsampled from a larger
#'   phylogenomic dataset based on their level of clock-likeness or phylogenetic
#'   signal, or otherwise at random. Each subsample was also run under two
#'   different models of molecular evolution, the site-homogeneous GTR+G and the
#'   site-heterogeneous CAT+GTR+G (Lartillot & Philippe 2004).
#'
#' @references
#' Mongiardino Koch, N., Thompson, J. R., Hiley, A. S., McCowin, M. F., Armstrong,
#'   A. F., Coppard, S. E., Aguilera, F., Bronstein, O., Kroh, A., Mooi, R., & Rouse,
#'   G. W. (2022). \emph{Phylogenomic analyses of echinoid diversification prompt a
#'   re-evaluation of their fossil record}. Elife, 11.
#' Lartillot N, Rodrigue N, Stubbs D, Richer J. 2013. PhyloBayes MPI: phylogenetic
#'   reconstruction with infinite mixtures of profiles in a parallel environment.
#'   Systematic Biology 62:611–615. DOI: https://doi.org/10.1093/sysbio/syt022,
#'   PMID: 23564032
"echinoid_dates"
