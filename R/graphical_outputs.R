#plot chronospace-------------------------------------------------------------------

#' Plot chronospace ordination(s) and chronogram warpings
#'
#' @importFrom ggtree %<+%
#' @import ggplot2
#'
#' @description For each factor, generate the two basic graphical
#'   representations of a chronospace: the projection of the sampled chronograms
#'   into the bgPCA ordination axes (a histogram if the factor has two levels, a
#'   scatter-plot otherwise), and the 'theoretical' representations of branch
#'   lengths for chronograms occupying the extremes of each axis.
#'
#' @param x An object of class \code{"chronospace"} containing one or more
#'   ordinations.
#' @param output Character, specifying which output should be plotted. Options
#'   are \code{"ordination"} for histograms / scatter plots of bgPC axes,
#'   \code{"extremes"} for trees showing node variation captured by bgPC axes,
#'   \code{"all"} for plotting both outputs. Option \code{"none"} is equivalent
#'   to \code{"all"}, yet allows results to be saved to an object without being
#'   plotted.
#' @param sdev Numeric, indicating at how many standard deviations should the
#'   extremes of the chronospace axes be depicted. Set to 1 sdev by default.
#'   Note that negative branch lengths might be generated, in which case the
#'   function automatically decreases the value until finding a representation
#'   for which all branch lengths are positive.
#' @param colors The colors used to represent groups (i.e. levels) of each
#'   factor.
#' @param factor Numeric; the factor or factors whose results are to be
#'   plotted and retained. If \code{NULL}, all the factors are retained.
#' @param axes Numeric of length two, specifying the axes to be depicted (not
#'   meaningful for factors with less than three levels).
#' @param ellipses Logical, indicating whether to plot data ellipses for each
#'   group (not meaningful for factors with less than three levels).
#' @param centroids Logical, indicating whether to plot group centroids (not
#'   meaningful for factors with less than three levels).
#' @param distances Logical, indicating whether to plot lines between group
#'   centroids. The width of the lines is inversely proportional to the
#'   distances between centroids in the original variable space (not meaningful
#'   for factors with less than three levels).
#' @param pt.alpha Numeric, indicating the transparency level of individual
#'   points in the scatter (not meaningful for factors with less than three
#'   levels).
#' @param pt.size Numeric, indicating the size of individual points in the
#'   scatter (not meaningful for factors with less than three levels).
#' @param ell.width Numeric, indicating line width for data ellipses (not
#'   meaningful for factors with less than three levels).
#' @param dist.width Numeric; scaling factor for the width of lines representing
#'   multivariate distances between group centroids (not meaningful for
#'   factors with less than three levels).
#' @param ct.size Numeric, indicating the size of the points marking group
#'   centroids (not meaningful for factors with less than three levels).
#' @param timemarks Numeric; an optional vector containing ages to be marked by
#'   vertical lines in chronospace representations.
#' @param gscale Logical; whether to add chronostratigraphic scale to trees
#'   (via \code{deeptime}).
#' @param ... Additional arguments passed to \code{plot}
#'
#' @details Starting from the object returned by [chronospace()], this function
#'   creates the two basic types of plots allowing interpretation of the
#'   ordination maximizing discrimination in node ages between groups (for
#'   each factor). The first of these is a projection of the chronograms
#'   (whose ages were used to generate the ordination) into the ordination axes,
#'   and can be either a histogram for two-level factors, or a bivariate
#'   scatter plot for factors with three or more levels. The second output
#'   contains trees representing the positive and negative extremes of each
#'   ordination axis, depicting the variation in nodes ages that it captures.
#'
#' @return A list containing the histogram/scatterplot and axes' extremes for
#'   each factor included.
#'
#' @export
#'
#' @references Mongiardino Koch N, Milla Carmona P (2024). Chronospaces: an R
#'   package for the statistical exploration of divergence times reveals extreme
#'   dependence on molecular clocks and gene choice. bioRxiv 2024.02.04.578835;
#'   doi: https://doi.org/10.1101/2024.02.04.578835.
#'
#' @examples
#' #Load ages data
#' data("echinoid_dates")
#'
#' #Create chronospace
#' cspace <- chronospace(echinoid_dates)
#'
#' #Plot chronospace ordination
#' csp.ord <- plot(cspace, output = "ordination")
#'
#' #Call same ordination for factor A from object
#' csp.ord$factor_A$ordination
#'
#' #Plot chronospace ordination
#' csp.ext <- plot(cspace, output = "extremes")
#'
#' #Call extremes of the first bgPC axis for factor A from object
#' csp.ext$factor_A$PC_extremes$bgPC1
#'
#' #Show univariate ordination for factor B from object
#' csp.ord$factor_B$ordination
#'
#' #Show extremes of the (only) bgPC axis for factor B object (notice all the
#' #results are stored there, even when output = "ordination")
#' csp.ord$factor_B$PC_extremes$bgPC1
plot.chronospace <- function(x, output = "all", sdev = 1, timemarks = NULL, gscale = TRUE,
                             ellipses = FALSE, centroids = FALSE, distances = FALSE,
                             colors = 1:5, factor = 1:(length(obj)), axes = c(1, 2),
                             pt.alpha = 0.5, pt.size = 1.5, ell.width = 1.2,
                             dist.width = 1, ct.size = 5, ...) {


  if(all(output != c("all", "ordination", "extremes", "none")))
    stop("output must be one of 'all', 'ordination', 'extremes', or 'none'")

  obj <- x[names(x) != "Total_vartable"]
  if(length(axes) != 2) axes <- c(1, 2)

  #create object for storing overall results, assign names
  results <- vector(mode = "list", length = length(factor))
  names(results) <- facnames <- names(obj)[factor]


  #get ordinations and PC extremes for factor i
  for(i in 1:length(factor)) {

    #create object for storing results of factor i, assign names
    results_i <- vector(mode = "list", length = 2)
    names(results_i) <- c("ordination", "PC_extremes")

    #extract information for factor i
    ordination <- obj[[factor[i]]]$ordination
    groups <- obj[[factor[i]]]$data$groups
    ages <- obj[[factor[i]]]$data$ages
    tree <- obj[[factor[i]]]$data$tree
    totvar <- obj[[factor[i]]]$ssq$totvar

    #set axes to either 1 (univariate plot) if the variable contains only two
    #groups, or 2 (bivariate plot) if it includes more groups
    num_functions <- 1
    if(ncol(ordination$x) >= 2) num_functions = 2

    #gather data for plotting
    colnames(ordination$x) <- NULL
    if(num_functions == 1) {
      to_plot <- data.frame(coordinates = ordination$x, groups = groups)
    } else {
      to_plot <- data.frame(coordinates = ordination$x[,axes], groups = groups)
    }

    #plot chronospace
    if(num_functions == 1) { #univariate

      chronospace <- ggplot(to_plot, aes(x = coordinates, fill = groups)) +
        geom_histogram(alpha = 0.5, position = 'identity', bins = 30) + theme_bw() +
        scale_fill_manual(values = colors) + ylab('Count') +
        theme(legend.title = element_blank(), panel.grid = element_blank()) +
        xlab(paste0('bgPCA axis 1 (',
                    round((100 * apply(ordination$x, 2, stats::var)[1] / totvar), 2), '% of variance)'))
    } else { #bivariate
      #compute groups centroids from bgPCA scores
      cents <- apply(X = ordination$x, MARGIN = 2, FUN = tapply, groups, mean)
      cents_df <- data.frame(coordinates.1 = cents[,axes[1]],
                             coordinates.2 = cents[,axes[2]],
                             groups = rownames(cents))

      #compute groups centroids from original variables;
      #calculate and standardize distances between centroids
      cents_original <- apply(X = ages, MARGIN = 2, FUN = tapply, groups, mean)
      dists <- as.matrix(stats::dist(cents_original))
      dists_std <- (1 / dists) / min(1 / dists)

      #generate combinations
      combins <- utils::combn(x = levels(groups), m = 2)

      #plot chronospace
      chronospace <- ggplot(to_plot, aes(x = coordinates.1, y = coordinates.2, color = groups)) +
        geom_point(alpha = pt.alpha, size = pt.size, key_glyph = "point") +
        theme_bw() + scale_color_manual(values = colors) +
        theme(legend.title = element_blank(), panel.grid = element_blank()) +
        xlab(paste0('bgPCA axis ', axes[1],  ' (',
                    round((100 * apply(ordination$x,2, stats::var)[axes[1]] / totvar), 2), '% of variance)')) +
        ylab(paste0('bgPCA axis ', axes[2],  ' (',
                    round((100 * apply(ordination$x,2, stats::var)[axes[2]] / totvar), 2), '% of variance)'))

      if(ellipses) {
        chronospace <- chronospace +
          stat_ellipse(lwd = ell.width, key_glyph = "point")
      }

      if(distances) {
        for(h in 1:ncol(combins)) {
          rdf <- cents_df[combins[,h],]
          width <- dists_std[combins[1,h], combins[2,h]]
          chronospace <- chronospace +
            geom_line(data = rdf, aes(x = coordinates.1, y = coordinates.2),
                      color = grDevices::gray.colors(n = 10)[1], linewidth = width * dist.width)
        }
      }

      if(centroids | distances) {
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

    #add title and save chronospace
    chronospace <- chronospace + ggtitle(paste0("Factor : ", facnames[i])) + theme(plot.title = element_text(hjust = 0.5))
    if(all(output != c("extremes", "none"))) print(chronospace)
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

      ax <- if(num_functions == 1) 1 else axes

      #loop through the bgPCA axes (depending on the number of groups in the
      #variable being tested)
      for(j in 1:num_functions) {

        #create a tree that contains topology but no branch lengths
        tree$edge.length <- rep(0, length(tree$edge.length))

        #constrain bgPC extremes to have positive branch lenghts only
        blen1 <- -1
        blen2 <- -1
        used_sdev <- sdev
        while(any(blen1 < 0) | any(blen2 < 0)) {

          #adjust bgPC range
          xrange <- c(used_sdev * stats::sd(tapply(ordination$x[,ax[j]], groups, mean)),
                      -used_sdev * stats::sd(tapply(ordination$x[,ax[j]], groups, mean)))

          #backwards PCA towards ages
          assign(paste0('plus_sd_', j), revPCA(xrange[1], ordination$rotation[,ax[j]], mean))
          assign(paste0('minus_sd_', j), revPCA(xrange[2], ordination$rotation[,ax[j]], mean))


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
            used_sdev <- round(used_sdev - 0.05, digits = 2)
          }
        }

        if(sdev != used_sdev) warning(paste0(facnames[i], " : sdev = ", sdev,
                                             " generated negative branch lengths for bgPC",
                                             ax[j], "; sdev = ", used_sdev, " was used instead."))

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

        if(all(output != c("ordination", "none"))) print(PCextremes[[j]])

      }

      #assign list names and save
      names(PCextremes) <- paste0("bgPC", axes[1:j])
      results_i$PC_extremes <- PCextremes

    }

    #add to overall results list
    results[[i]] <- results_i

  }

  return(invisible(results))

}


#get user-specified nodes ------------------------------------------------------

#' Plot the posterior distribution of ages of a user-specified node
#'
#' @import ggplot2
#'
#' @description Plot the distribution of a given node across inference
#'   conditions.
#'
#' @param data_ages An object of class \code{"nodeAges"}.
#' @param tips Character vector of length 2, specifying the tip names of two
#'   terminals bracketing the node whose age (expressed in million of years) is
#'   to be plotted.
#' @param factor Numeric; the factor or factors whose results are to be
#'   plotted and retained. If \code{NULL}, all the factors are retained.
#' @param plot whether to plot the results, or only store them (default=TRUE)
#' @param colors The colors used to represent groups (i.e. levels) of each
#'   factor.
#' @param l.alpha Numeric, indicating the transparency level of density lines.
#' @param l.width Numeric, indicating the width (thickness) of density lines.
#' @param timemarks Numeric; an optional vector containing ages to be marked by
#'   vertical lines.
#' @param gscale Logical; whether to add chronostratigraphic scale to trees
#'   (via \code{deeptime}).
#'
#' @details This function takes a single character vector containing the tip
#'   names of two terminals, identifies their most recent common ancestor, and
#'   plots the distribution of posterior ages for said node across the
#'   conditions explored.
#'
#' @return A panel showing the distribution of ages for the target node under
#'   each level of the requested factors.
#'
#' @export
#'
#' @references
#'
#' @examples
#' #Load ages data
#' data("echinoid_dates")
#'
#' #Get ages distribution for the MRCA of Brissus obesus and Abatus cordatus
#' MRCA_Brissus_Abatus <- specified_node(echinoid_dates, tips = c('Brissus_obesus',
#'    'Abatus_cordatus'), plot = FALSE)
#'
#' #Show age distribution for the MRCA of the two terminals associated with
#' #factor A
#' MRCA_Brissus_Abatus$factor_A
specified_node <- function(data_ages, tips = NULL, factor = 1:ncol(data_ages$factors),
                           plot = TRUE, colors = 1:5, l.alpha = 0.3, l.width = 2,
                           timemarks = NULL, gscale = FALSE) {

  #create data_ages object for storing overall results, assign names
  results <- vector(mode = "list", length = length(factor))
  names(results) <- colnames(data_ages$factors)[factor]

  #for each factor:
  for(i in 1:length(factor)) {

    #extract information for factor i
    groups <- data_ages$factors[,factor[i]]
    ages <- data_ages$ages
    tree <- data_ages$topology

    #find ages corresponding to node of interest
    #to make sure this is consistent with how the nodeAges object was built,
    #we'll use the same approach
    clades <- list()
    for(j in 1:tree$Nnode) {
      clades[j] <- list(tree$tip.label[unlist(phangorn::Descendants(tree, length(tree$tip.label) + j, type = 'tips'))])
    }

    #which of these nodes include the two species
    includes_targets <- c()
    for(j in 1:length(clades)) {
      if(all(tips %in% clades[[j]])) {
        includes_targets <- c(includes_targets, j)
      }
    }

    #the least inclusive node that includes both targets is their MRCA
    target_clade_size <- min(sapply(clades, length)[includes_targets])
    node_number <- includes_targets[which(sapply(clades, length)[includes_targets] ==
                                           target_clade_size)]

    #restrict ages to just that one node
    ages <- ages[,colnames(ages) == paste0('clade_', node_number)]
    gmeans <- apply(as.data.frame(ages), 2, tapply, INDEX = groups, mean)

    #plot the posterior distributions of the chosen node

    #make room to save the individual plots
    plots <- vector(mode = "list", length = 1)

    #plot
    timemarks1.1 <- timemarks[timemarks <= max(ages) & timemarks >= min(ages)]
    to_plot <- data.frame(age = ages, group = groups)
    plots <- ggplot(to_plot, aes(x = age, color = group)) +
      geom_density(alpha = l.alpha, linewidth = l.width) + scale_x_reverse() +
      theme_bw() + scale_color_manual(values = colors) +
      theme(plot.title = element_text(size = 8), panel.grid = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
      xlab('Age of MRCA') + ylab('Density') +
      geom_vline(xintercept = timemarks1.1, lty = 2, col = "gray")

    if(gscale) {
    plots <- plots +
      deeptime::coord_geo(size = 3.5, height = unit(1, "line"),
                          expand = TRUE, pos = "bottom")
    }

    if(length(unique(to_plot$group)) == 2) {
      plots <- plots +
        ggtitle(paste0('MRCA of ', tips[1], ' and ',
                       tips[2], ' (difference = ',
                       round((max(gmeans) - min(gmeans)), 1),
                       ' Ma)'))
    } else {
      plots <- plots +
        ggtitle(paste0('MRCA of ', tips[1], ' and ',
                       tips[2], ' (max difference = ',
                       round((max(gmeans) - min(gmeans)), 1),
                       ' Ma)'))
    }

    #plot and save
    if(plot) print(plots)

    results[[i]] <- plots

  }

  return(invisible(results))
}


#get sensitive nodes ----------------------------------------------------

#' Identify the most sensitive nodes and depict their age distribution
#'
#' @import ggplot2
#'
#' @description For each factor, identify the most variables nodes (in terms of
#'   age) for each factor, and plot their distribution of posterior ages.
#'
#' @param data_ages An object of class \code{"nodeAges"}.
#' @param amount_of_change Numeric, specifying the desired amount of variation
#'   in node age (expressed in million of years) above which the nodes are
#'   retained and depicted.
#' @param num_clades Numeric, indicating the desired number of most sensitive
#'   nodes to be retained and depicted. Ignored if \code{amount_of_change} is
#'   specified.
#' @param factor Numeric; the factor or factors whose results are to be
#'   plotted and retained. If \code{NULL}, all the factors are retained.
#' @param plot whether to plot the results, or only store them (default=TRUE)
#' @param colors The colors used to represent groups (i.e. levels) of each
#'   factor.
#' @param l.alpha Numeric, indicating the transparency level of density lines.
#' @param l.width Numeric, indicating the width (thickness) of density lines.
#' @param timemarks Numeric; an optional vector containing ages to be marked by
#'   vertical lines.
#' @param gscale Logical; whether to add chronostratigraphic scale to trees
#'   (via \code{deeptime}).
#'
#' @details This function identifies, for each factor, the nodes in the fixed
#'   topology whose ages are most sensitive (i.e. variable). This can be done
#'   either by indicating a threshold for variation in age, above which nodes
#'   are considered 'sensitive', or by specifying the number of most sensitive
#'   nodes which must be retained. The nodes are named using two terminals
#'   selected at random from the two descendant subclades of each node.
#'
#'   For each of these most sensitive nodes, the function will plot the
#'   posterior distribution of ages by level.
#'
#' @return A panel showing the posterior distribution for the age of
#'   each of the most sensitive nodes associated with each factor.
#'
#' @export
#'
#' @examples
#' #Load ages data
#' data("echinoid_dates")
#'
#' #Get the 5 most sensitive nodes
#' sensinodes5 <- sensitive_nodes(echinoid_dates, num_clades = 5, plot = FALSE)
#'
#' #Show ages distribution for the 5 most sensitive nodes associated to factor A
#' sensinodes5$factor_A
sensitive_nodes <- function(data_ages, amount_of_change = NULL, num_clades = 5,
                            factor = 1:ncol(data_ages$factors), plot = TRUE,
                            colors = 1:5, l.alpha = 0.3, l.width = 2,
                            timemarks = NULL, gscale = FALSE) {

  #create data_ages object for storing overall results, assign names
  results <- vector(mode = "list", length = length(factor))
  names(results) <- colnames(data_ages$factors)[factor]

  #for each factor:
  for(i in 1:length(factor)) {

    #extract information for factor i
    groups <- data_ages$factors[,factor[i]]
    ages <- data_ages$ages
    tree <- data_ages$topology

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
      if(is.null(num_clades)) {
        num_nodes <- 5
      } else {
        #else go with what the user chose, although cap at 20
        num_nodes <- num_clades
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
      this_ages <- ages[,which(colnames(ages) == names(clade))]
      timemarks1.1 <- timemarks[timemarks <= max(this_ages) & timemarks >= min(this_ages)]
      to_plot <- data.frame(age = this_ages, group = groups)
      plots[[j]] <- ggplot(to_plot, aes(x = age, color = group)) +
        geom_density(alpha = l.alpha, linewidth = l.width) + scale_x_reverse() +
        theme_bw() + scale_color_manual(values = colors) +
        theme(plot.title = element_text(size = 8), panel.grid = element_blank(),
              axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
        xlab('Age of MRCA') + ylab('Density') +
        geom_vline(xintercept = timemarks1.1, lty = 2, col = "gray")

      if(gscale) {
        plots[[j]] <- plots[[j]] +
          deeptime::coord_geo(size = 3.5, height = unit(1, "line"),
                              expand = TRUE, pos = "bottom")
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
    if(num_clades <= 5) {
      row_number = num_clades
    } else {
      row_number = 5
    }
    most_affected <- ggpubr::annotate_figure(ggpubr::ggarrange(plotlist = plots,
                                                               common.legend = TRUE, legend = 'bottom',
                                                               ncol = ceiling(num_nodes / 5), nrow = row_number))
    if(plot) print(most_affected)

    results[[i]] <- most_affected
  }

  return(invisible(results))
}


#LTT by group-------------------------------------------------------------------

#' Plot average lineage-through-time (LTT) curves
#'
#' @importFrom magrittr %>%
#' @import ggplot2
#'
#' @description For each factor, plot the LTT curve obtained under each level by
#'   averaging across the corresponding chronograms.
#'
#' @param data_ages A \code{"nodeAges"} object created using [extract_ages()].
#' @param summary Character, indicating whether the 'mean' or 'median' is to be
#'   used in computations.
#' @param uncertainty Character, indicating whether (and how to) represent age
#'   uncertainty. Options are 'none' (default), 'CI_<number>' (with <number>
#'   being a value between 1 and 99, representing the extent of a confidence
#'   interval, e.g. 'CI_95'), or 'sample_<number>' (with <number> being a number
#'   of random posterior typologies to be drawn, e.g. 'sample_20').
#' @param colors The colors used to represent groups (i.e. levels) of each
#'   factor.
#' @param l.alpha Numeric, indicating the transparency level of density lines.
#' @param l.width Numeric, indicating the width (thickness) of density lines.
#' @param factor Numeric; the factor or factors whose results are to be
#'   plotted and retained. If \code{NULL}, all the factors are retained.
#' @param plot whether to plot the results, or only store them (default=TRUE)
#' @param timemarks Numeric; an optional vector containing ages to be marked by
#'   vertical lines in LTT plots.
#' @param gscale Logical; whether to add chronostratigraphic scale to trees
#'   (via \code{deeptime}).
#'
#' @export
#'
#' @references Harvey PH, May RM, & Nee S. (1994). Phylogenies without fossils.
#'   Evolution, 48: 523–529.
#'
#' @examples
#' #Load ages data
#' data("echinoid_dates")
#'
#' #Create LTT plots
#' sensi_ltt <- ltt_sensitivity(echinoid_dates, summary = "mean", plot = FALSE)
#'
#' #Show LTT plot for factor A only
#' sensi_ltt$factor_A
ltt_sensitivity <- function(data_ages, summary = 'median', uncertainty = 'none',
                            colors = 1:5, l.alpha = 0.5, l.width = 2,
                            factor = 1:ncol(data_ages$factors),
                            plot = TRUE, timemarks = NULL, gscale = FALSE) {

  ages <- data_ages$ages
  factors <- data_ages$factors

  ltts <- vector(mode = "list", length = length(factor))
  names(ltts) <- colnames(factors[factor])

  for(i in 1:length(factor)) {

    this_groups <- factors[,factor[i]]

    if(summary == 'mean') {
      ages_average <- ages %>% dplyr::mutate(type = this_groups) %>%
        dplyr::group_by(type) %>%
        dplyr::summarise_if(is.numeric, mean, .groups = 'drop') %>%
        dplyr::mutate(metric = 'average') %>%
        dplyr::select(type, metric, tidyr::everything())
    }

    if(summary == 'median') {
      ages_average <- ages %>% dplyr::mutate(type = this_groups) %>%
        dplyr::group_by(type) %>%
        dplyr::summarise_if(is.numeric, stats::median, .groups = 'drop') %>%
        dplyr::mutate(metric = 'average') %>%
        dplyr::select(type, metric, tidyr::everything())
    }

    if(grepl('CI', uncertainty)) {
      level <- as.numeric(gsub('CI_', '', uncertainty)) / 100
      level <- (1 - level) / 2

      for(j in 1:length(levels(this_groups))) {
        this_ages <- ages[this_groups == levels(this_groups)[j],]
        this_ages <- apply(this_ages, 2, function(x) sort(x, decreasing  = F))

        uncertain_ages <- cbind(type = rep(levels(this_groups)[j], 2),
                                metric = c('min', 'max'),
                                data.frame(this_ages[c(ceiling((nrow(this_ages) * level) + 1),
                                            floor((nrow(this_ages) * (1 - level)) - 1)),]))

        if(j == 1) {
          all_uncertain <- uncertain_ages
        } else {
          all_uncertain <- rbind(all_uncertain, uncertain_ages)
        }

      }
    }

    if(grepl('sample', uncertainty)) {
      sample <- as.numeric(gsub('sample_', '', uncertainty))

      for(j in 1:length(levels(this_groups))) {
        this_ages <- ages[this_groups == levels(this_groups)[j],]
        this_ages <- this_ages[sample(1:nrow(this_ages), sample, replace = F),]

        uncertain_ages <- cbind(type = rep(levels(this_groups)[j], sample),
                                metric = paste0(levels(this_groups)[j], '_sample_', 1:sample),
                                data.frame(this_ages))

        if(j == 1) {
          all_uncertain <- uncertain_ages
        } else {
          all_uncertain <- rbind(all_uncertain, uncertain_ages)
        }
      }
    }

    if(exists('all_uncertain')) {
      ages_average <- rbind(ages_average, all_uncertain)
    }

    ages_average <- dplyr::arrange(ages_average, type)

    ages_average[,3:ncol(ages_average)] <- t(apply(ages_average[,3:ncol(ages_average)],
                                                   1,
                                                   sort, decreasing = T))

    ages_average <- tidyr::pivot_longer(ages_average,
                                        cols = starts_with('clade'),
                                        names_to = 'num_lineages',
                                        values_to = 'age') %>%
      dplyr::mutate(num_lineages = as.numeric(gsub('clade_', '',
                                                   num_lineages)) + 1)

    to_add <- ages_average %>% dplyr::mutate(num_lineages = num_lineages - 1)
    ages_average <- rbind(ages_average, to_add,
                          expand.grid(type = levels(ages_average$type),
                                      metric = unique(ages_average$metric),
                                      num_lineages = max(ages_average$num_lineages),
                                      age = 0)) %>%
      dplyr::arrange(type, num_lineages)

    timemarks1.1 <- timemarks[timemarks <= max(ages_average$age) &
                                timemarks >= min(ages_average$age)]

    ltts[[i]] <- ggplot(ages_average, aes(color = type)) +
      scale_color_manual(values = colors) +
      geom_line(data = subset(ages_average, ages_average$metric == 'average'),
                aes(x = age, y = num_lineages, color = type),
                alpha = l.alpha, linewidth = l.width) +
      scale_y_log10() + scale_x_reverse() +
      theme_bw() + xlab('Age (Ma)') + ylab('Number of lineages') +
      theme(panel.grid = element_blank()) +
      geom_vline(xintercept = timemarks1.1, lty = 2, col = "gray") +
      guides(colour = guide_legend(override.aes = list(alpha = 1, shape = 21,
                                                       fill = colors[1:nlevels(this_groups)],
                                                       size = 3.5)))

    if(grepl('CI', uncertainty)) {
      ltts[[i]] <- ltts[[i]] +
        scale_fill_manual(values = colors) +
        geom_ribbon(data = tidyr::unnest(tidyr::pivot_wider(subset(ages_average,
                                                                   ages_average$metric != 'average'),
                                                            names_from = 'metric',
                                                            values_from = 'age', values_fn = list),
                                         cols = c(min, max)) %>%
                      dplyr::arrange(type, num_lineages),
                    aes(xmin = min, xmax = max, y = num_lineages, fill = type),
                    color = NA, alpha = 0.1)
    }

    if(grepl('sample', uncertainty)) {
      ltts[[i]] <- ltts[[i]] +
        geom_line(data = subset(ages_average,
                                ages_average$metric != 'average') %>%
                    dplyr::arrange(type, metric, num_lineages),
                  aes(x = age, y = num_lineages, color = type, group = metric),
                  alpha = 0.1)
    }

    if(gscale) {
      ltts[[i]] <- ltts[[i]] +
        deeptime::coord_geo(size = 3.5, height = unit(1, "line"), expand = TRUE, pos = "bottom")
    }

    if(plot) print(ltts[[i]])

  }

  return(invisible(ltts))
}
