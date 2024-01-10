#plot chronospace-------------------------------------------------------------------

#' Plot chronospace ordination(s) and axes' extremes
#'
#' @importFrom ggtree %<+%
#' @import ggplot2
#'
#' @description For each factor, generate the two basic graphical
#'   representations of a chronospace: the projection of the sampled chronograms
#'   into the ordination axes, and the two 'theoretical' extremes of those axes.
#'
#' @param obj An object of class \code{"chronospace"} containing one or more
#'   ordinations.
#' @param sdev Numeric, indicating at how many standard deviations should the
#'   extremes of the chronospace axes be depicted. If NULL (the default),
#'   extremes will be depicted at the highest standard deviation avoiding
#'   negative branch lengths.
#' @param colors The colors used to represent groups (i.e. levels) of each
#'   factor.
#' @param factors Numeric; the factor or factors whose results are to be
#'   retained (by default, the first two axes).
#' @param axes Numeric of length two, specifying the axes to be depicted (not
#'   meaningful for factors with less than three levels).
#' @param ellipses Logical, indicating whether to plot data ellipses for each
#'   group (not meaningful for factors with less than three levels).
#' @param centroids Logical, indicating whether to plot groups' centroids (not
#'   meaningful for factors with less than three levels).
#' @param distances Logical, indicating whether to plot lines between groups'
#'   centroids whose width is proportional to the distances between them in the
#'   original variable space (not meaningful for factors with less than three
#'   levels).
#' @param pt.alpha Numeric, indicating the transparency level of individual
#'   points in the scatter (not meaningful for factors with less than three
#'   levels).
#' @param pt.size Numeric, indicating the size of individual points in the
#'   scatter (not meaningful for factors with less than three levels).
#' @param ell.width Numeric, indicating line width for data ellipses (not
#'   meaningful for factors with less than three levels).
#' @param dist.width Numeric; scaling factor for the width of lines representing
#'   multivariate distances between groups' centroids (not meaningful for
#'   factors with less than three levels).
#' @param ct.size Numeric, indicating the size of the points marking groups'
#'   centroids (not meaningful for factors with less than three levels).
#' @param timemarks Numeric; an optional vector containing ages to be marked by
#'   vertical lines in chronospace representations.
#' @param gscale Logical; whether to add chronostratigraphic scale to trees
#'   (via \code{deeptime}).
#'
#' @details Starting from the object returned by [chronospace()], this function
#'   creates the two basic types of plots allowing interpretation of the
#'   ordination maximizing variation in node age between the groups (of
#'   each factor). The first of these is a projection of the phylogenetic trees
#'   (whose ages were used to generate the ordination) into the ordination axes,
#'   and can be either a histogram for two-level factors, or a bivariate
#'   scatterplot for factors with three or more levels. The second output
#'   contains 'theoretical' trees representing the positive and negative
#'   extremes of each ordination axis, depicting the variation in nodes ages
#'   it captures.
#'
#' @return A list containing the histogram/scatterplot and axes' extremes for
#'   each factor included.
#'
#' @export
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
#' #Plot chronospace
#' plot(cspace)
#'
#' #Show (bivariate) ordination for factor A
#' cspace$factor_A$ordination
#'
#' #Show extremes of the first bgPC axis for factor A
#' cspace$factor_A$PC_extremes
#'
#' #' #Show (univariate) ordination for factor B
#' cspace$factor_B$ordination
#'
#' #Show extremes of the bgPC axis for factor B
#' cspace$factor_B$PC_extremes
plot.chronospace <- function(obj, sdev = NULL, timemarks = NULL, gscale = TRUE,
                             ellipses=TRUE, centroids=FALSE, distances = FALSE,
                             colors = 1:5, factors = 1:length(obj), axes = c(1, 2),
                             pt.alpha = 0.5, pt.size = 1.5, ell.width = 1.2,
                             dist.width = 1, ct.size = 5) {


  #obj <- obj[names(obj) != "Total"] #########################
  obj <- obj[names(obj) != "Total_vartable"]

  if(length(axes) != 2) axes <- c(1, 2)

  #create object for storing overall results, assign names
  results <- vector(mode = "list", length = length(obj))
  names(results) <- facnames <- names(obj)

  warns <- NULL

  #get ordinations and PC extremes for factor i
  for(i in 1:length(obj)){

    #create object for storing results of factor i, assign names
    results_i <- vector(mode = "list", length = 2)
    names(results_i) <- c("ordination", "PC_extremes")

    #extract information for factor i
    #######################################
    # bgPCA <- obj[[i]]
    # groups <- bgPCA$groups
    # totvar <- bgPCA$totvar
    # ages <- bgPCA$ages
    # tree <- bgPCA$tree
    #######################################
    ordination <- obj[[i]]$ordination
    groups <- obj[[i]]$data$groups
    ages <- obj[[i]]$data$ages
    tree <- obj[[i]]$data$tree
    totvar <- obj[[i]]$ssq$totvar

    #set axes to either 1 (univariate plot) if the variable contains only two
    #groups, or 2 (bivariate plot) if it includes more groups
    num_functions <- 1
    #if(ncol(bgPCA$x) >= 2) num_functions = 2  ##############
    if(ncol(ordination$x) >= 2) num_functions = 2

    #gather data for plotting
    #colnames(bgPCA$x) <- NULL ##########
    colnames(ordination$x) <- NULL
    if(num_functions == 1) {
      #to_plot <- data.frame(coordinates = bgPCA$x, groups = groups) ############
      to_plot <- data.frame(coordinates = ordination$x, groups = groups)
    } else {
      #to_plot <- data.frame(coordinates = bgPCA$x[,axes], groups = groups) ##########
      to_plot <- data.frame(coordinates = ordination$x[,axes], groups = groups)
    }

    #plot chronospace
    if(num_functions == 1) { #univariate

      chronospace <- ggplot(to_plot, aes(x = coordinates, fill = groups)) +
        geom_histogram(alpha = 0.5, position = 'identity', bins = 30) + theme_bw() +
        scale_fill_manual(values = colors) + ylab('Count') +
        theme(legend.title = element_blank(), panel.grid = element_blank()) +
        xlab(paste0('bgPCA axis 1 (',
                    #round((100 * apply(bgPCA$x, 2, stats::var)[1] / totvar), 2), '% of variance)')) #########
                    round((100 * apply(ordination$x, 2, stats::var)[1] / totvar), 2), '% of variance)'))
    } else { #bivariate
      #compute groups centroids from bgPCA scores
      #cents <- apply(X = bgPCA$x, MARGIN = 2, FUN = tapply, groups, mean) #########
      cents <- apply(X = ordination$x, MARGIN = 2, FUN = tapply, groups, mean)
      cents_df <- data.frame(coordinates.1 = cents[,axes[1]],
                             coordinates.2 = cents[,axes[2]],
                             groups = rownames(cents))

      #compute groups centroids from original variables;
      #calculate and standardize distances between centroids
      cents_original <- apply(X = ages, MARGIN = 2, FUN = tapply, groups, mean)
      dists <- as.matrix(stats::dist(cents_original))
      dists_std <- dists / max(dists)

      #generate combinations
      combins <- utils::combn(x = levels(groups), m = 2)

      #plot chronospace
      chronospace <- ggplot(to_plot, aes(x = coordinates.1, y = coordinates.2, color = groups)) +
        geom_point(alpha = pt.alpha, size = pt.size, key_glyph = "point") +
        theme_bw() + scale_color_manual(values = colors) +
        theme(legend.title = element_blank(), panel.grid = element_blank()) +
        xlab(paste0('bgPCA axis ', axes[1],  ' (',
                    #round((100 * apply(bgPCA$x,2, stats::var)[axes[1]] / totvar), 2), '% of variance)')) + ########
                    round((100 * apply(ordination$x,2, stats::var)[axes[1]] / totvar), 2), '% of variance)')) +
        ylab(paste0('bgPCA axis ', axes[2],  ' (',
                    #round((100 * apply(bgPCA$x,2, stats::var)[axes[2]] / totvar), 2), '% of variance)')) #######
                    round((100 * apply(ordination$x,2, stats::var)[axes[2]] / totvar), 2), '% of variance)'))

      if(ellipses){
        chronospace <- chronospace +
          stat_ellipse(lwd = ell.width, key_glyph = "point")
      }

      if(distances){
        for(h in 1:ncol(combins)){
          rdf <- cents_df[combins[,h],]
          width <- (5 * dists_std[combins[1,h], combins[2,h]]) - 2
          chronospace <- chronospace +
            geom_line(data = rdf, aes(x = coordinates.1, y = coordinates.2),
                      color = grDevices::gray.colors(n = 10)[1], size = width * dist.width)
        }
      }

      if(centroids | distances){
        chronospace <- chronospace +
          geom_point(data = cents_df, aes(x = coordinates.1, y = coordinates.2),
                     color = "black", shape = 21,
                     fill = colors[1:nlevels(groups)],
                     size = ct.size)
      }

      chronospace <- chronospace +
        guides(colour = guide_legend(override.aes = list(alpha = 1, shape = 21,
                                                         color = "black",
                                                         fill = colors[1:nlevels(groups)],
                                                         size = 3.5)))
    }

    #save chronospace
    print(chronospace + ggtitle(paste0("Factor : ", facnames[i])) + theme(plot.title = element_text(hjust = 0.5)))
    results_i$ordination <- chronospace

    #Finally, compute changes in each branch captured by the bgPCA axes (needs a
    #tree!)
    if(is.na(tree)[1]) {
      cat('Plotting changes on branch lengths can only be shown if a tree is provided\n')
    } else {
      #obtain clades from tree
      clades <- list()
      for(j in 1:tree$Nnode) {
        clades[j] <- list(tree$tip.label[unlist(phangorn::Descendants(tree, length(tree$tip.label) + j,
                                                                      type = 'tips'))])
      }

      #ages implied by a position at the origin of the bgPCA plot
      mean <- matrix(colMeans(ages), ncol = 1)

      #create object for storing the extremes of the bgPC j
      PCextremes <- vector(mode = "list", length = num_functions)

      if(num_functions == 1) {
        ax <- 1
      } else {
        ax <- axes
      }

      #loop through the bgPCA axes (depending on the number of groups in the
      #variable being tested)
      for(j in 1:num_functions) {

        #create a tree that contains topology but no branch lengths
        tree$edge.length <- rep(0, length(tree$edge.length))

        #if sdev has been specified, use them to obtain extremes; otherwise, constrain
        #bgPC extremes to have positive branch lenghts only
        if(!is.null(sdev)) {
          ######################################################################################
          # xrange <- c(sdev * stats::sd(tapply(bgPCA$x[,ax[j]], groups, mean)),
          #             -sdev * stats::sd(tapply(bgPCA$x[,ax[j]], groups, mean)))
          # assign(paste0('plus_sd_', j),  revPCA(xrange[1], bgPCA$rotation[,ax[j]], mean))
          # assign(paste0('minus_sd_', j), revPCA(xrange[2], bgPCA$rotation[,ax[j]], mean))
          ######################################################################################

          xrange <- c(sdev * stats::sd(tapply(ordination$x[,ax[j]], groups, mean)),
                      -sdev * stats::sd(tapply(ordination$x[,ax[j]], groups, mean)))
          assign(paste0('plus_sd_', j),  revPCA(xrange[1], ordination$rotation[,ax[j]], mean))
          assign(paste0('minus_sd_', j), revPCA(xrange[2], ordination$rotation[,ax[j]], mean))


          #retrieve trees with branch lenghts corresponding to ages at the extremes of bgPC j
          extrees <- reconstruct_blen(clades = clades,
                                      tree = tree,
                                      plus = get(paste0('plus_sd_', j)),
                                      mean = mean,
                                      minus = get(paste0('minus_sd_', j)))

          tree_plus <- extrees$tree_plus
          tree_mean <- extrees$tree_mean
          tree_minus <- extrees$tree_minus

          used_sdev <- sdev

        } else {
          adj <- 1
          blen1 <- -1
          blen2 <- -1
          while(any(blen1 < 0) | any(blen2 < 0)) {

            #adjust bgPC range
            ########################################################################
            # xrange <- c(stats::sd(tapply(bgPCA$x[,ax[j]], groups, mean)),
            #             -stats::sd(tapply(bgPCA$x[,ax[j]], groups, mean)))
            ########################################################################
            xrange <- c(stats::sd(tapply(ordination$x[,ax[j]], groups, mean)),
                        -stats::sd(tapply(ordination$x[,ax[j]], groups, mean)))

            newmax <- xrange[1] + (diff(xrange) * (1 - adj) / 2)
            newmin <- xrange[2] - (diff(xrange) * (1 - adj) / 2)

            #backwards PCA towards ages
            ########################################################################
            # assign(paste0('plus_sd_', j), revPCA(newmax, bgPCA$rotation[,ax[j]], mean))
            # assign(paste0('minus_sd_', j), revPCA(newmin, bgPCA$rotation[,ax[j]], mean))
            ########################################################################
            assign(paste0('plus_sd_', j), revPCA(newmax, ordination$rotation[,ax[j]], mean))
            assign(paste0('minus_sd_', j), revPCA(newmin, ordination$rotation[,ax[j]], mean))

            #get % of sd used to get only positive branch lenghts
            # used_sdev <- round(newmax / stats::sd(tapply(bgPCA$x[,ax[j]], groups, mean)), 3) ######
            used_sdev <- round(newmax / stats::sd(tapply(ordination$x[,ax[j]], groups, mean)), 3)

            #retrieve trees with branch lenghts corresponding to ages at the extremes of bgPC j
            extrees <- reconstruct_blen(clades = clades,
                                        tree = tree,
                                        plus = get(paste0('plus_sd_', j)),
                                        minus = get(paste0('minus_sd_', j)),
                                        mean = mean)

            tree_plus <- extrees$tree_plus
            tree_mean <- extrees$tree_mean
            tree_minus <- extrees$tree_minus

            blen1 <- tree_plus$edge.length
            blen2 <- tree_minus$edge.length
            if(any(blen1 < 0) | any(blen2 < 0)) {
              adj <- adj - 0.05
            }
          }
        }

        #compute delta in branch lengths between the mean tree and the positive and negative extremes
        changes_plus <- tree_plus$edge.length - tree_mean$edge.length
        changes_minus <- tree_minus$edge.length - tree_mean$edge.length


        #if time marks have been specified, use them to  draw vertical lines in the corresponding tree
        if(!is.null(timemarks)){
          t.range <- range(phytools::nodeHeights(tree_minus))
          timemarks1.1 <- timemarks[timemarks <= t.range[2] & timemarks >= t.range[1]]

          t.range <- range(phytools::nodeHeights(tree_plus))
          timemarks2.1 <- timemarks[timemarks <= t.range[2] & timemarks >= t.range[1]]
        } else {
          timemarks1.1 <- timemarks2.1 <- NULL
        }


        #convert phylo trees into ggtrees, adding delta in branch length to the metadata
        tree_plus_gg <- suppressMessages(ggtree::ggtree(tree_plus, size = 1.5) %<+%
                                           data.frame(node = tree_plus$edge[,2], delta = changes_plus))
        tree_minus_gg <- suppressMessages(ggtree::ggtree(tree_minus, size = 1.5) %<+%
                                            data.frame(node = tree_minus$edge[,2], delta = changes_minus))

        warn <- if(any(na.omit(tree_minus_gg$data$branch.length / abs(tree_minus_gg$data$branch.length)) == -1)) TRUE else FALSE
        warns <- c(warns, warn)



        #create graphics for each extreme of the bgPC j
        tree_minus_gg$data$x <- max(tree_minus_gg$data$x) - tree_minus_gg$data$x
        negative <- tree_minus_gg + aes(color=delta) +
          scale_color_gradient2(limits = range(c(changes_minus, changes_plus)),
                                high = "red", low = "blue", mid = "gray", midpoint = 0) +
          ggtitle(paste0(facnames[i], " - bgPC", axes[j], ", \nnegative extreme (sd = ", -used_sdev,")")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_vline(xintercept = timemarks1.1, lty = 2, col = "gray")

        tree_plus_gg$data$x <- max(tree_plus_gg$data$x) - tree_plus_gg$data$x
        positive <- tree_plus_gg + aes(color = delta) +
          scale_color_gradient2(limits = range(c(changes_minus, changes_plus)),
                                high = "red", low = "blue", mid = "gray", midpoint = 0) +
          ggtitle(paste0(facnames[i], " - bgPC", axes[j], ", \npositive extreme (sd = ", used_sdev,")")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_vline(xintercept = timemarks2.1, lty = 2, col = "gray")

        #add chronostratigraphic scale
        if(gscale == TRUE) {
          negative <- negative +
            scale_x_reverse() +
            deeptime::coord_geo(size = 3.5, height = unit(1, "line"), expand = TRUE)

          positive <- positive +
            scale_x_reverse() +
            deeptime::coord_geo(size = 3.5, height = unit(1, "line"), expand = TRUE)
        }

        #combine both into a single graphic and store
        PCextremes[[j]] <- negative + positive + patchwork::plot_layout(guides = "collect") &
          theme(legend.position = "bottom")

      }

      #assign list names and save
      names(PCextremes) <- paste0("bgPC", axes[1:j])
      results_i$PC_extremes <- PCextremes

    }

    #add to overall results list
    results[[i]] <- results_i

  }

  if(any(warns)) warning(paste("sdev =", sdev, "generates negative branch lengths"))

  factors <- factors[!factors > length(results)]
   return(results[factors])

}


#get sensitive nodes ----------------------------------------------------


#' Identify the most sensitive nodes and depict their age distribution
#'
#' @import ggplot2
#'
#' @description For each node, identify the most variables nodes (in terms of
#'   age), and plot their age proportion distributions.
#'
#' @param obj An object of class \code{"nodeAges"}.
#' @param amount_of_change Numeric, specifying the desired amount of variation
#'   in node age (expressed in million of years) above which the nodes are
#'   retained and depicted.
#' @param chosen_clades Numeric, indicating the desired number of most sensitive
#'   nodes to be retained and depicted.
#' @param factors Numeric; the factor or factors whose results are to be
#'   retained (by default, the first two axes). Ignored if amount_of_change is
#'   specified.
#' @param colors The colors used to represent groups (i.e. levels) of each
#'   factor.
#' @param timemarks Numeric; an optional vector containing ages to be marked by
#'   vertical lines.
#' @param gscale Logical; whether to add chronostratigraphic scale to trees
#'   (via \code{deeptime}).
#'
#' @details This function identifies, for each factor, the nodes in the fixed
#'   topology whose ages are most sensitive (i.e. variable) as a function of the
#'   corresponding levels. These can be done either by indicating a threshold
#'   for variation in age above which nodes are considered relevant, or by
#'   specifying the number of most sensitive nodes which must be retained. The
#'   nodes are named using two terminals selected at random from the two
#'   subclades defined by each node.
#'
#'   For each of these most sensitive nodes, the function will plot the
#'   distribution of relative proportion of ages, by level.
#'
#' @return A panel showing the relative proportions distribution for the age of
#'   each of the most sensitive nodes associated to each factor.
#'
#' @export
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
#' #Get the 5 most sensitive nodes
#' sensinodes5 <- sensitive_nodes(obj = cspace, chosen_clades = 5)
#'
#' #Show ages distribution for the 5 most sensitive nodes associated to factor A
#' sensinodes5$factor_A
sensitive_nodes <- function(obj, amount_of_change = NULL, chosen_clades = 5, factors = 1:ncol(obj$factors),
                            colors = 1:5, timemarks = NULL, gscale = FALSE) {


  #create object for storing overall results, assign names
  results <- vector(mode = "list", length = ncol(obj$factors))
  names(results) <- colnames(obj$factors)


  #for each factor:
  for(i in 1:ncol(obj$factors)) {

    #extract information for factor i
    groups <- obj$factors[,i]
    ages <- obj$ages
    tree <- obj$topology

    #compute mean age of each level
    gmeans <- apply(ages, 2, tapply, INDEX = groups, mean)

    #plot the posterior distribution of nodes with the strongest differences
    #between runs

    #first decide how many nodes will be plotted
    #if an minimum amount of change is specified, go with it
    if(!is.null(amount_of_change)) {
      num_nodes <- length(which((apply(gmeans, 2, max) -
                                   apply(gmeans, 2, min)) > amount_of_change))

      #reduce to a max of 20,or plot 5 if none changes by the specified amount
      if(num_nodes > 20) num_nodes <- 20
      if(num_nodes == 0) num_nodes <- 5

    } else { #if a minimum amount is not specified
      #if a number of clades is not specified, do 5
      if(is.null(chosen_clades)) {
        num_nodes <- 5
      } else {
        #else go with what the user chose, although cap at 20
        num_nodes <- chosen_clades
      }
      if(num_nodes > 20) num_nodes <- 20
    }

    #obtain clades from tree
    clades = list()
    for(j in 1:tree$Nnode) {
      clades[j] <- list(tree$tip.label[unlist(phangorn::Descendants(tree, length(tree$tip.label) + j,
                                                                    type = 'tips'))])
    }

    #make room to save the individual plots
    plots <- vector(mode = "list", length = num_nodes)

    #loop through the nodes
    for(j in 1:num_nodes) {

      #sort clades starting by those that vary the most between analyses and
      #choose clade j
      clade <- which(sort((apply(gmeans, 2, max) - apply(gmeans, 2, min)), decreasing = TRUE)[j] ==
                       (apply(gmeans, 2, max) - apply(gmeans, 2, min)))


      #obtain corresponding node number and the descendant taxa
      node <- phangorn::mrca.phylo(tree, clades[[clade]])
      desc <- phangorn::Descendants(tree, node = node, type = 'children')

      #obtain representative taxa from either side of the split
      for(k in 1:length(desc)) {
        #if it is a tip, extract the name
        if(as.numeric(desc[k]) <= length(tree$tip.label)) {
          desc[k] <- tree$tip.label[as.numeric(desc[k])]
          #else choose a random tip from the descendant clade
        } else {
          desc[k] <- tree$tip.label[sample(unlist(phangorn::Descendants(tree, node = as.numeric(desc[k]),
                                                                        type = 'tips')), 1)]
        }
      }

      #plot
      ages_clade_scaled <- (max(ages[,clade]) - ages[,clade]) + min(ages[,clade])
      timemarks1.1 <- timemarks[timemarks <= max(ages_clade_scaled) & timemarks >= min(ages_clade_scaled)]
      to_plot <- data.frame(age = ages_clade_scaled, group = groups)
      plots[[j]] <- ggplot(to_plot, aes(x = age, color = group)) +
        geom_density(alpha = 0.3, size = 2) +
        theme_bw() + scale_color_manual(values = colors) +
        theme(plot.title = element_text(size = 8), panel.grid = element_blank()) +
        scale_x_continuous(breaks = pretty(to_plot$age), labels = abs(pretty(to_plot$age))) +
        xlab('Age of MRCA') + ylab('Density') +
        geom_vline(xintercept = timemarks1.1, lty = 2, col = "gray")

      if(gscale == TRUE) {
        plots[[j]] <- plots[[j]] +
          deeptime::coord_geo(size = 3.5, height = unit(1, "line"), expand = TRUE, pos = "bottom")
      }

      if(length(unique(to_plot$group)) == 2) {
        plots[[j]] <- plots[[j]] +
          ggtitle(paste0('MRCA of ', desc[1], ' and ',
                         desc[2], ' (difference = ',
                         round((max(gmeans[,clade]) - min(gmeans[,clade])), 1),
                         ' Ma)'))
      } else {
        plots[[j]] <- plots[[j]] +
          ggtitle(paste0('MRCA of ', desc[1], ' and ',
                         desc[2], ' (max difference = ',
                         round((max(gmeans[,clade]) - min(gmeans[,clade])), 1),
                         ' Ma)'))
      }
    }

    #plot and save, accounting for a varying number of columns depending on the
    #nodes plotted
    most_affected <- ggpubr::annotate_figure(ggpubr::ggarrange(plotlist = plots,
                                                               common.legend = T, legend = 'bottom',
                                                               ncol = ceiling(num_nodes / 5), nrow = 5))
    results[[i]] <- most_affected
  }

  factors <- factors[!factors > length(results)]
  return(results[factors])
}



#LTT by group-------------------------------------------------------------------

#' Plot average Lineage Through Time (LTT) curves
#'
#' @importFrom magrittr %>%
#' @import ggplot2
#'
#' @description For each factor, plot the LTT curve of each level averaged
#'   across the corresponding subsample of chronograms.
#'
#' @param data_ages A \code{"nodeAges"} object created using [extract_ages()].
#' @param average Character, indicating whether the 'mean' or 'median' is to be
#'   used in computations.
#' @param colors The colors used to represent groups (i.e. levels) of each
#'   factor.
#' @param timemarks Numeric; an optional vector containing ages to be marked by
#'   vertical lines in LTT plots.
#' @param gscale Logical; whether to add chronostratigraphic scale to trees
#'   (via \code{deeptime}).
#' @export
#'
#' @examples
#' #Load ages data
#' data("data_ages")
#'
#' #Create LTT plots
#' sensiltt <- ltt_sensitivity(data_ages = data_ages, average = "mean")
#'
#' #Show LTT plot for factor A
#' sensiltt$factor_A
ltt_sensitivity <- function(data_ages, average = 'median', colors = 1:5,
                            timemarks = NULL, gscale = TRUE) {

  ages <- data_ages$ages
  factors <- data_ages$factors

  ltts <- vector(mode = "list", length = ncol(factors))
  names(ltts)<-colnames(factors)

  for(i in 1:ncol(factors)) {

    sample <- nrow(factors)/length(unique(factors[,i]))
    num_nodes <- ncol(ages)

    this_ages <- apply(ages, 1, sort)
    this_groups <- factors[,i]
    this_order <- order(this_groups)

    this_ages <- this_ages[,this_order]
    this_groups <- this_groups[this_order]

    colnames(this_ages) <- 1:ncol(this_ages)
    this_ages <- tidyr::pivot_longer(tibble::as_tibble(this_ages), 1:ncol(this_ages)) %>%
      dplyr::mutate(name = as.numeric(name)) %>% dplyr::arrange(name, dplyr::desc(value)) %>%
      dplyr::mutate(type = rep(as.character(unique(this_groups)),
                               each = sample * num_nodes),
                    num_lineages = rep(2:(num_nodes + 1), length(this_groups)))

    if(average == 'mean') {
      ages_average <- this_ages %>% dplyr::group_by(type, num_lineages) %>%
        dplyr::summarise(av_value = mean(value), .groups = 'drop')
    }
    if(average == 'median') {
      ages_average <- this_ages %>% dplyr::group_by(type, num_lineages) %>%
        dplyr::summarise(av_value = stats::median(value), .groups = 'drop')
    }

    to_add <- ages_average %>% dplyr::mutate(num_lineages = num_lineages - 1)
    ages_average <- rbind(ages_average, to_add) %>% dplyr::arrange(type, num_lineages)
    timemarks1.1 <- timemarks[timemarks <= max(ages_average$av_value) &
                                timemarks >= min(ages_average$av_value)]

    ltts[[i]] <- ggplot(ages_average, aes(x = av_value, y = num_lineages, color = type)) +
      scale_color_manual(values = colors) +
      geom_line(alpha = 0.5, size = 2) + scale_y_log10() + scale_x_reverse() +
      theme_bw() + xlab('Age (Ma)') + ylab('Number of lineages') +
      theme(panel.grid = element_blank()) +
      geom_vline(xintercept = timemarks1.1, lty = 2, col = "gray") +
      guides(colour = guide_legend(override.aes = list(alpha = 1, shape = 21,
                                                       fill = colors[1:nlevels(this_groups)],
                                                       size = 3.5)))

    if(gscale == TRUE) {
      ltts[[i]] <- ltts[[i]] +
        deeptime::coord_geo(size = 3.5, height = unit(1, "line"), expand = TRUE, pos = "bottom")
    }

  }

  return(ltts)

}

