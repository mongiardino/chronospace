#create chronospace-------------------------------------------------------------------

#' Create chronospace
#'
#' @description Compute the ordination maximizing variation in node ages using
#'   between-group PCA (one for each factor).
#'
#' @param data_ages A matrix created using [extract_ages()].
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
chronospace <- function(data_ages, vartype = "non-redundant")  {

  #split data.frame 'data_ages' into ages and factors
  ages <- data_ages[,which(grepl('clade', colnames(data_ages)))]
  groups <- data_ages[,which(grepl('factor', colnames(data_ages)))]

  #factors names
  if(is.null(dim(groups))) groups <- data.frame(groups)
  facnames <- paste0("factor_", LETTERS[1:ncol(groups)])

  #create object for storing overall results, assign names
  results <- vector(mode = "list", length = ncol(groups))
  names(results) <- facnames

  #perform bgPCA using each factor separately
  for(i in 1:ncol(groups)) {

    #perform bgPCA between groups defined by factor i over original variation
    bgPCA1 <- bgprcomp(x = ages, groups = groups[,i])

    #compute percentage of variation explained
    totvar <- sum(apply(ages, 2, var))
    expvar <- sum(apply(bgPCA1$x, 2, var))
    perc_tot <- 100 * (expvar / totvar)

    #report proportion of total variation explained
    if(ncol(groups) > 1) {
      catnames <- paste0('(', paste(unique(groups[,i]), collapse = '/'), ')')
      cat(paste('--- Results for', facnames[i], catnames, '---\n'))
    } else {
      cat('Results:\n')
    }

    cat(paste0('Proportion of total variation in node ages explained by ',
               facnames[i], ' = ',
               round(perc_tot, digits = 3),
               '%', '\n'))

    if(ncol(groups) > 1){
      #use bgPCA to compute an ordination that is residual to all factors but factor i
      bgPCA2.1 <- bgprcomp(x = ages, groups = groups[,-i])
      resids2.1 <- bgPCA2.1$residuals

      #perform bgPCA between groups defined by factor i over residual variation
      bgPCA2.2 <- bgprcomp(x = resids2.1, groups = groups[,i])
      expvar2.2 <- sum(apply(bgPCA2.2$x, 2, var))

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
    if(vartype == "total" | ncol(groups) == 1) bgPCA <- bgPCA1
    if(vartype == "non-redundant" & ncol(groups) > 1) bgPCA <- bgPCA2.2

    #store bgPCA results, along with total variation and groups of factor i
    bgPCA$totvar <- totvar
    bgPCA$groups <- groups[,i]
    bgPCA$ages <- ages
    results[[i]] <- bgPCA
  }

  class(results) <- "chronospace"

  return(invisible(results))
}

#create chronospace-------------------------------------------------------------------

#' Create chronospace
#'
#' @description Compute the ordination maximizing variation in node ages using
#'   between-group PCA (one for each factor).
#'
#' @param data_ages A matrix created using [extract_ages()].
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
chronospace2 <- function(data_ages, vartype = "non-redundant")  {

  #split data.frame 'data_ages' into ages and factors
  ages <- data_ages$ages
  factors <- data_ages$factors

  #factors names
  if(is.null(dim(factors))) factors <- data.frame(factors)
  facnames <- paste0("factor_", LETTERS[1:ncol(factors)])

  #create object for storing overall results, assign names
  results <- vector(mode = "list", length = ncol(factors))
  names(results) <- facnames

  #perform bgPCA using each factor separately
  for(i in 1:ncol(factors)) {

    #perform bgPCA between groups defined by factor i over original variation
    bgPCA1 <- bgprcomp(x = ages, groups = factors[,i])

    #compute percentage of variation explained
    totvar <- sum(apply(ages, 2, var))
    expvar <- sum(apply(bgPCA1$x, 2, var))
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
      bgPCA2.1 <- bgprcomp(x = ages, groups = factors[,-i])
      resids2.1 <- bgPCA2.1$residuals

      #perform bgPCA between groups defined by factor i over residual variation
      bgPCA2.2 <- bgprcomp(x = resids2.1, groups = factors[,i])
      expvar2.2 <- sum(apply(bgPCA2.2$x, 2, var))

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
    if(vartype == "total" | ncol(factors) == 1) bgPCA <- bgPCA1
    if(vartype == "non-redundant" & ncol(factors) > 1) bgPCA <- bgPCA2.2

    #store bgPCA results, along with total variation and groups of factor i
    bgPCA$totvar <- totvar
    bgPCA$groups <- factors[,i]
    bgPCA$ages <- ages
    bgPCA$tree <- data_ages$topology
    results[[i]] <- bgPCA
  }

  class(results) <- "chronospace"

  return(invisible(results))
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

