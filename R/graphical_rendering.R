#plot chronospace-------------------------------------------------------------------

#' Plot chronospace ordination(s) and axes' extremes
#'
#' @description For each factor, generate the two basic graphical representations of a chronospace: the projection of
#' the sampled chronograms into the synthetic chronospace axes, and the two ‘theoretical’ extremes of those axes.
#'
#' @param obj An object containing one or more ordinations created using [chronospace()].
#' @param tree An object of class "phylo" containing the same fixed topology as the trees sampled from the posterior.
#' @param sdev Numeric, indicating at how many standard deviations should the extremes of the chronospace axes be
#' depicted.
#' @param timemarks Numeric; an optional vector containing ages to be marked by vertical lines in chronospace
#' representations.
#' @param colors The colors used to represent groups (i.e. levels) of each factor.
#' @param factors Numeric; the factor or factors whose results are to be retained (by default, the first two axes).
#' @param axes Numeric of length two, specifying the axes to be depicted (not meaningful for factors with less than
#' three levels).
#' @param ellipses Logical, indicating whether to plot data ellipses for each group (not meaningful for factors with
#' less than three levels).
#' @param centroids Logical, indicating whether to plot groups' centroids (not meaningful for factors with less than
#' three levels).
#' @param distances Logical, indicating whether to plot lines between groups' centroids whose width is proportional
#' to the distances between them in the original variable space (not meaningful for factors with less than three levels).
#' @param pt.alpha Numeric, indicating the transparency level of individual points in the scatter (not meaningful for
#' factors with less than three levels).
#' @param pt.size Numeric, indicating the size of individual points in the scatter (not meaningful for factors with
#' less than three levels).
#' @param ell.width Numeric, indicating line width for data ellipses (not meaningful for factors with less than three
#' levels).
#' @param dist.width Numeric; scaling factor for the width of lines representing multivariate distances between groups'
#' centroids (not meaningful for factors with less than three levels).
#' @param ct.size Numeric, indicating the size of the points marking groups' centroids (not meaningful for factors with
#' less than three levels).
#'
#' @details Starting from the object returned by [chronospace()], this function creates the two basic types of plots
#' allowing interpretation of the synthetic ordination maximizing variation in node age between the groups of each factor.
#' The first of these is a graphical projection of the samples of phylogenetic trees used to generate the ordination into
#' its synthetic axes, which can be either a histogram for factors with only two levels, or a bivariate scatterplot for
#' factors with three or more levels. The second output consists of ‘theoretical’ trees representing the positive and
#' negative extremes of each synthetic axis, depicting the variation in nodes ages captured by it.
#'
#' @return A list containing the histogram/scatterplot and axes' extremes for each factor included.
#' @export
#'
#' @examples
plot.chronospace<-function(obj, tree=NA, sdev=1, timemarks = NULL,
                           ellipses=TRUE, centroids=FALSE, distances=FALSE,
                           colors=1:5, factors=1:length(obj), axes=c(1,2), pt.alpha=0.5, pt.size=1.5, ell.width=1.2, dist.width=1, ct.size=5) {

  if(length(axes)!=2) axes<-c(1,2)

  #create object for storing overall results, assign names
  results <- vector(mode = "list", length = length(obj))
  names(results) <- facnames <- names(obj)

  #get ordinations and Pc extremes for factor i
  for(i in 1:length(obj)){

    #create object for storing results of factor i, assing names
    results_i <- vector(mode = "list", length = 2)
    names(results_i) <- c("ordination", "PC_extremes")

    #extract information for factor i
    bgPCA <- obj[[i]]
    groups <- bgPCA$groups
    totvar <- bgPCA$totvar
    ages <- bgPCA$ages

    #set axes to either 1 (univariate plot) if the variable contains only two
    #groups, or 2 (bivariate plot) if it includes more groups
    num_functions <- 1
    if(ncol(bgPCA$x) >= 2) num_functions = 2

    #gather data for plotting
    to_plot <- data.frame(coordinates = bgPCA$x, groups = groups)

    #plot chronospace
    if(num_functions == 1) { #univariate

      chronospace <- ggplot(to_plot, aes(x = coordinates, fill = groups)) +
        geom_histogram(alpha = 0.5, position = 'identity', bins = 30) + theme_bw() +
        scale_fill_manual(values = colors) + ylab('Count') +
        theme(legend.title = element_blank(), panel.grid = element_blank()) +
        xlab(paste0('bgPCA axis 1 (', round((100*apply(bgPCA$x,2,var)[1]/totvar), 2), '% of variance)'))

    } else { #bivariate
      #compute groups centroids from bgPCA scores
      cents<-apply(X = bgPCA$x, MARGIN = 2, FUN = tapply, groups, mean)
      cents_df<-data.frame(coordinates.1=cents[,1], coordinates.2=cents[,2], groups=rownames(cents))

      #compute groups centroids from original variables; calculate and standardize distances between centroids
      cents_original<-apply(X = ages, MARGIN = 2, FUN = tapply, groups, mean)
      dists<-as.matrix(dist(cents_original))
      dists_std<-dists/max(dists)

      #generate combinations
      combins<-combn(x = levels(groups), m = 2)

      #plot chronospace
      chronospace<-ggplot(to_plot, aes(x = coordinates.1, y = coordinates.2, color = groups)) +
        geom_point(alpha = pt.alpha, size=pt.size, key_glyph = "point") +
        theme_bw() + scale_color_manual(values = colors) +
        theme(legend.title = element_blank(), panel.grid = element_blank()) +
        xlab(paste0('bgPCA axis ', axes[1],  ' (', round((100*apply(bgPCA$x,2,var)[1]/totvar), 2), '% of variance)')) +
        ylab(paste0('bgPCA axis ', axes[2],  ' (', round((100*apply(bgPCA$x,2,var)[2]/totvar), 2), '% of variance)'))

      if(ellipses){
        chronospace <- chronospace +
          stat_ellipse(lwd=ell.width, key_glyph = "point")
      }

      if(distances){
        for(h in 1:ncol(combins)){
          rdf<-cents_df[combins[,h],]
          width<-(5*dists_std[combins[1,h], combins[2,h]])-2
          chronospace <- chronospace + geom_line(data=rdf, aes(x = coordinates.1, y = coordinates.2), color=gray.colors(n=10)[1], size=width*dist.width)
        }
      }

      if(centroids|distances){
        chronospace <- chronospace +
          geom_point(shape=21, data=cents_df, color="black", fill = colors[1:nlevels(groups)], aes(x = coordinates.1, y = coordinates.2), size=ct.size)
      }

      chronospace <- chronospace +
        guides(colour = guide_legend(override.aes = list(alpha=1, shape=21, color="black", fill = colors[1:nlevels(groups)], size=3.5)))

    }

    #save chronospace
    results_i$ordination <- chronospace

    #obtain clades from tree
    clades <- list()
    for(j in 1:tree$Nnode) {
      clades[j] <- list(tree$tip.label[unlist(Descendants(tree, length(tree$tip.label)+j, type = 'tips'))])
    }

    #Finally, compute changes in each branch captured by the bgPCA axes (needs a
    #tree!)
    if(is.na(tree)[1]) {
      cat('Plotting changes on branch lengths can only be shown if a tree is provided\n')
    } else {
      #ages implied by a position at the origin of the bgPCA plot
      mean <- matrix(colMeans(ages), ncol = 1)

      #create object for storing the extremes of the bgPC j
      PCextremes <- vector(mode = "list", length = num_functions)

      if(num_functions==1) ax<-1 else ax<-axes

      #loop through the bgPCA axes (depending on the number of groups in the
      #variable being tested)
      for(j in 1:num_functions) {

        #create a tree that contains topology but no branch lengths
        tree$edge.length <- rep(0, length(tree$edge.length))

        #ages implied by moving along this bgPCA axis 'sdev' number of standard
        #deviations to both sides
        assign(paste0('plus_sd_', j), revPCA(sdev*sd(bgPCA$x[,ax[j]]), bgPCA$rotation[,ax[j]], mean))
        assign(paste0('minus_sd_', j), revPCA(-sdev*sd(bgPCA$x[,ax[j]]), bgPCA$rotation[,ax[j]], mean))

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
              node_to_change <- getMRCA(tree, unlist(clades[which_clades[l]]))

              #get node ages for this node
              dif_minus <- get(paste0('minus_sd_', j))[,which_clades[l]]
              dif_mean <- mean[which_clades[l],]
              dif_plus <- get(paste0('plus_sd_', j))[,which_clades[l]]

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
                nodes_of_descendants <- tree$edge[,2][which(tree$edge[,1] == node_to_change)]

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
                  tips <- unlist(Descendants(tree, nodes_of_descendants[m],
                                             type = 'tips'))

                  #remove from the age the age of the descendant node, which is
                  #already set up correctly as the loop goes from smaller to
                  #larger clades
                  tree_minus$edge.length[which(tree_minus$edge[,2] == nodes_of_descendants[m])] <-
                    dif_minus - dist.nodes(tree_minus)[tips[1], nodes_of_descendants[m]]
                  tree_mean$edge.length[which(tree_mean$edge[,2] == nodes_of_descendants[m])] <-
                    dif_mean - dist.nodes(tree_mean)[tips[1], nodes_of_descendants[m]]
                  tree_plus$edge.length[which(tree_plus$edge[,2] == nodes_of_descendants[m])] <-
                    dif_plus - dist.nodes(tree_plus)[tips[1], nodes_of_descendants[m]]
                }
              }
            }
          }
        }

        #compute delta in branch lengths between the mean tree and the positive and negative extremes
        changes_plus <- tree_plus$edge.length - tree_mean$edge.length
        changes_minus <- tree_minus$edge.length - tree_mean$edge.length

        #if time marks have been specified, use them to  draw vertical lines in the corresponding tree
        if(!is.null(timemarks)){
          t.max <- max(nodeHeights(tree_minus))
          timemarks1.1 <- timemarks[which(timemarks <= t.max)]
          timemarks1.2 <- t.max - timemarks1.1

          t.max <- max(nodeHeights(tree_plus))
          timemarks2.1 <- timemarks[which(timemarks <= t.max)]
          timemarks2.2 <- t.max-timemarks2.1
        } else {
          timemarks1.2 <- timemarks2.2 <- NULL
        }

        #convert phylo trees into ggtrees, adding delta in branch length to the metadata
        tree_plus_gg <- ggtree(tree_plus, size = 1.5) %<+%
          data.frame(node = tree_plus$edge[,2], delta = changes_plus)
        tree_minus_gg <- ggtree(tree_minus, size = 1.5) %<+%
          data.frame(node = tree_minus$edge[,2], delta = changes_minus)

        #create graphics for each extreme of the bgPC j
        negative <- tree_minus_gg + aes(color=delta) +
          scale_color_gradient2(limits = range(c(changes_minus, changes_plus)),
                                high = "red", low = "blue", mid = "gray", midpoint = 0) +
          ggtitle(paste0(facnames[i], " - bgPC", j, ", negative extreme")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_vline(xintercept = timemarks1.2, lty = 2, col = "gray")

        positive <- tree_plus_gg + aes(color = delta) +
          scale_color_gradient2(limits = range(c(changes_minus, changes_plus)),
                                high = "red", low = "blue", mid = "gray", midpoint = 0) +
          ggtitle(paste0(facnames[i], " - bgPC", j, ", positive extreme")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_vline(xintercept = timemarks2.2, lty = 2, col = "gray")

        #combine both into a single graphic and store
        PCextremes[[j]] <- negative + positive + plot_layout(guides = "collect") &
          theme(legend.position="bottom")

      }

      #assign list names and save
      names(PCextremes) <- paste0("bgPC", 1:j)
      results_i$PC_extremes <- PCextremes

    }

    #add to overall results list
    results[[i]] <- results_i

  }

  factors<-factors[!factors>length(results)]
  return(results[factors])

}


#get senstive nodes ----------------------------------------------------

#' Obtain the most sensitive nodes and depict their age distribution
#'
#' @description Identify the most sensitive nodes associated to each factor, and plot their ages proportion distributions.
#'
#' @param obj An object containing one or more ordinations created using [chronospace()].
#' @param tree An object of class "phylo" containing the same fixed topology as the trees from the posterior.
#' @param amount_of_change Numeric, specyfing the desired amount of variation in age (expressed in million of years)
#' above which the nodes are retained and depicted.
#' @param chosen_clades Numeric, indicating the desired number of most sensitive nodes to be retained and depicted.
#' @param factors Numeric; the factor or factors whose results are to be retained (by default, the first two axes).
#' Ignored if amount_of_change is explicited.
#' @param colors The colors used to represent groups (i.e. levels) of each factor.
#'
#' @details This function identifies, for each factor, the nodes in the fixed topology whose ages are most sensitive
#' (i.e. variable) as a function of the corresponding levels. These can be done either by indicating a threshold for
#' variation in age above which nodes are considered relevant, or by specifying the number of most sensitive nodes which
#' must be retained. The nodes are named using two terminals selected at random from the two subclades defined by each
#' node.
#'
#' For each of these most sensitive nodes, the function will plot the distribution of relative proportion of ages, by level.
#'
#' @return A panel showing the relative proportions distribution for the age of each of the most sensitive nodes
#' #associated to each factor.
#' @export
#'
#' @examples
sensitive_nodes <- function(obj, tree, amount_of_change, chosen_clades, factors=1:length(obj), colors=1:5){

  #create object for storing overall results, assign names
  results <- vector(mode = "list", length = length(obj))
  names(results) <- names(obj)

  #perform bgPCA on each variable
  for(i in 1:length(obj)) {

    #extract information for factor i
    bgPCA <- obj[[i]]
    groups <- bgPCA$groups
    ages <- bgPCA$ages

    #plot the posterior distribution of nodes with the strongest differences
    #between runs

    #first decide how many nodes will be plotted
    #if an minimum amount of change is specified, go with it
    if(!is.na(amount_of_change)) {
      num_nodes <- length(which((apply(bgPCA$gmeans, 2, max) -
                                   apply(bgPCA$gmeans, 2, min)) > amount_of_change))

      #reduce to a max of 20,or plot 5 if none changes by the specified amount
      if(num_nodes > 20) num_nodes <- 20
      if(num_nodes == 0) num_nodes <- 5

    } else { #if a minimum amount is not specified
      #if a number of clades is not specified, do 5
      if(is.na(chosen_clades)) {
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
      clades[j] <- list(tree$tip.label[unlist(Descendants(tree, length(tree$tip.label)+j,
                                                          type = 'tips'))])
    }

    #make room to save the individual plots
    plots <- vector(mode = "list", length = num_nodes)

    #loop through the nodes
    for(j in 1:num_nodes) {

      #sort clades starting by those that vary the most between analyses and
      #choose clade j
      clade <- which(sort((apply(bgPCA$gmeans, 2, max) - apply(bgPCA$gmeans, 2, min)), decreasing = T)[j] ==
                       (apply(bgPCA$gmeans, 2, max) - apply(bgPCA$gmeans, 2, min)))


      #obtain corresponding node number and the descendant taxa
      node <- mrca.phylo(tree, clades[[clade]])
      desc <- Descendants(tree, node = node, type = 'children')

      #obtain representative taxa from either side of the split
      for(k in 1:length(desc)) {
        #if it is a tip, extract the name
        if(as.numeric(desc[k]) <= length(tree$tip.label)) {
          desc[k] <- tree$tip.label[as.numeric(desc[k])]
          #else choose a random tip from the descendant clade
        } else {
          desc[k] <- tree$tip.label[sample(unlist(Descendants(tree, node = as.numeric(desc[k]),
                                                              type = 'tips')), 1)]
        }
      }

      #make the plot
      to_plot <- data.frame(age = ages[,clade], group = groups)
      plots[[j]] <- ggplot(to_plot, aes(x = -age, color = group)) +
        geom_density(alpha = 0.3, size = 2) +
        theme_bw() + scale_color_manual(values = colors) +
        theme(plot.title = element_text(size = 8)) +
        scale_x_continuous(breaks = pretty(-to_plot$age), labels = abs(pretty(-to_plot$age))) +
        xlab('Age of MRCA') + ylab('Density')

      if(length(unique(to_plot$group)) == 2) {
        plots[[j]] <- plots[[j]] +
          ggtitle(paste0('MRCA of ', desc[1], ' and ',
                         desc[2], ' (difference = ',
                         round((max(bgPCA$gmeans[,clade]) - min(bgPCA$gmeans[,clade])), 1),
                         ' Ma)'))
      } else {
        plots[[j]] <- plots[[j]] +
          ggtitle(paste0('MRCA of ', desc[1], ' and ',
                         desc[2], ' (max difference = ',
                         round((max(bgPCA$gmeans[,clade])-min(bgPCA$gmeans[,clade])), 1),
                         ' Ma)'))
      }
    }

    #plot and save, accounting for a varying number of columns depending on the
    #nodes plotted
    most_affected <- annotate_figure(ggarrange(plotlist = plots,
                                               common.legend = T, legend = 'bottom',
                                               ncol = ceiling(num_nodes/5), nrow = 5))
    #plot(most_affected)
    results[[i]] <- most_affected
  }

  factors<-factors[!factors>length(results)]
  return(results[factors])
}


#LTT by group-------------------------------------------------------------------

#' Plot average Lineage Through Time (LTT) curves
#'
#' @description For each factor, plot the LTT curve of each level averaged across the corresponding  subsample of
#' chronograms.
#'
#' @param data_ages A matrix created using [extract_ages()].
#' @param average Character, indicating whether the 'mean' or 'median' is to be used in computations.
#'
#' @export
#'
#' @examples
ltt_sensitivity <- function(data_ages, average = 'median') {
  ages <- data_ages[,which(grepl('clade', colnames(data_ages)))]
  groups <- data_ages[,which(grepl('factor', colnames(data_ages)))]
  plots <- vector(mode = "list", length = ncol(groups))

  for(i in 1:ncol(groups)) {
    sample <- nrow(groups)/length(unique(groups[,i]))
    num_nodes <- ncol(ages)

    this_ages <- apply(ages, 1, sort)
    this_groups <- groups[,i]
    this_order <- order(this_groups)

    this_ages <- this_ages[,this_order]
    this_groups <- this_groups[this_order]

    colnames(this_ages) <- 1:ncol(this_ages)
    this_ages <- pivot_longer(as.tibble(this_ages), 1:ncol(this_ages)) %>%
      mutate(name = as.numeric(name)) %>% arrange(name, desc(value)) %>%
      mutate(type = rep(as.character(unique(this_groups)),
                        each = sample * num_nodes),
             num_lineages = rep(2:(num_nodes + 1), length(this_groups)))

    if(average == 'mean') {
      ages_average <- this_ages %>% group_by(type, num_lineages) %>%
        summarise(av_value = mean(value), .groups = 'drop')
    }
    if(average == 'median') {
      ages_average <- this_ages %>% group_by(type, num_lineages) %>%
        summarise(av_value = median(value), .groups = 'drop')
    }

    to_add <- ages_average %>% mutate(num_lineages = num_lineages - 1)
    ages_average <- rbind(ages_average, to_add) %>% arrange(type, num_lineages)

    plots[[i]] <- ggplot(ages_average, aes(x = av_value, y = num_lineages, color = type)) +
      geom_line(alpha = 0.3, size = 2) + scale_y_log10() + scale_x_reverse() +
      theme_bw() + xlab('Age (Ma)') + ylab('Number of lineages')
  }

  ltts <- annotate_figure(ggarrange(plotlist = plots,
                                    common.legend = F, legend = 'bottom',
                                    ncol = ncol(groups), nrow = 1))

  return(ltts)

}
