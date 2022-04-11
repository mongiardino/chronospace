# internal between-group PCA function ---------------------------------------------
bgprcomp <- function(x, groups){

  grandmean <- colMeans(x)
  x_centered <- scale(x, scale = F, center = T)
  x_gmeans <- apply(X = x_centered, MARGIN = 2, FUN = tapply, groups, mean)

  V_g <- cov(x_gmeans)
  eig <- eigen(V_g)

  scores <- x_centered%*%eig$vectors
  scores <- cbind(scores[,1:(nlevels(groups) - 1)])
  rotation <- eig$vectors

  preds <- scores %*% t(rotation[,1:ncol(scores)])
  resids <- x - preds

  return(list(x = scores, residuals = resids, rotation = rotation,
              values = eig$values, center = grandmean, gmeans = x_gmeans))
}


# internal reverse PCA function -----------------------------------------------------
revPCA<-function(scores, vectors, center){ t(t(scores%*%t(vectors))+center) }


#create chronospace-------------------------------------------------------------------
chronospace <- function(data_ages, variation = "non-redundant")  {

  #split data.frame 'data_ages' into ages and factors
  ages <- data_ages[,which(grepl('clade', colnames(data_ages)))]
  groups <- data_ages[,which(grepl('factor', colnames(data_ages)))]

  #factors names
  if(is.null(dim(groups))) groups <- data.frame(groups)
  facnames<-paste0("factor_", LETTERS[1:ncol(groups)])

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
    perc_tot <- 100 * (expvar/totvar)

    #report proportion of total variation explained
    cat(paste0('Proportion of total variation in node ages explained by ',
               facnames[i], ' = ',
               round(perc_tot, digits=3),
               '%', '\n'))

    if(ncol(groups)>1){
      #use bgPCA to compute an ordinaion that is residual to all factors but factor i
      bgPCA2.1 <- bgprcomp(x = ages, groups = groups[,-i])
      resids2.1 <- bgPCA2.1$residuals

      #perform bgPCA between groups defined by factor i over residual variation
      bgPCA2.2 <- bgprcomp(x = resids2.1, groups = groups[,i])
      expvar2.2 <- sum(apply(bgPCA2.2$x, 2, var))

      #compute percentage of non-redundant variation explained
      perc_nonred <- 100 * (expvar2.2/totvar)

      #report proportion of non-redundant variation explained
      cat(paste0('Proportion of non-redundant variation in node ages explained by ',
                 facnames[i], ' = ',
                 round(perc_nonred, digits=3),
                 '%', '\n'))
    } else {cat('(There is only one factor, non-redundant variation omitted)\n')}

    #select which bgPCA results are going to be used
    if(variation == "total" | ncol(groups)==1) bgPCA <- bgPCA1
    if(variation == "non-redundant" & ncol(groups)>1) bgPCA <- bgPCA2.2

    #store bgPCA results, along with total variation and groups of factor i
    bgPCA$totvar<-totvar
    bgPCA$groups<-groups[,i]
    bgPCA$ages<-ages
    results[[i]]<-bgPCA
  }

  return(invisible(results))
}
