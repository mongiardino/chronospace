
## `chronospace`: statistical exploration of time-calibrated phylogenies and the relative impact of methodological decisions.

<img src="man/figures/chronospace_hex.png" align="right" width="200"/>

Calibrating phylogenies to absolute time is a complex and
time-consuming, yet crucial step in macroevolutionary analysis. Bayesian
approaches to time scaling phylogenetic trees have dramatically grown in
complexity, and hinge today upon numerous methodological choices (e.g.,
gene sampling strategy, type of clock, model of molecular evolution,
prior on divergence times, set of fossil calibrations, etc.) that can be
difficult to justify. As a consequence, divergence times are routinely
inferred under a limited range of parametric conditions, often chosen
arbitrarily. Attempts to measure and summarize the sensitivity of
results to such decisions often rely on tables with node ages or
plotting multiple consensus trees which are cumbersome to interpret, and
lack a consistent way to quantify (and therefore to compare) their
effect.

`chronospace` is an R package devised to help visualizing and
quantifying the sensitivity of results to such decisions through the use
of chronospaces, i.e., graphical representations summarizing variation
in the node ages contained in time-calibrated trees with fixed topology.
In particular, `chronospace` uses between-groups PCA (bgPCA) for
summarizing variation in node ages produced by specific methodological
choices, and estimates their impact using a Sum of squares (SSQ)
approach for measuring effect size.

**Note: this repository addresses the version of `chronospace` described
in Mongiardino Koch and Milla Carmona (in preparation). If you are
interested in the data and procedures described in Mongiardino Koch et
al. 2022 eLife paper (on which this version of the package is based), go
to the \[chronospaces_eLife repository\]
(<https://github.com/mongiardino/chronospaces_eLife/tree/main>)**

## Installation

The development version of `chronospace` can be installed from
[GitHub](https://github.com/) using:

``` r
# install.packages("devtools")
devtools::install_github("mongiardino/chronospace")
```

## The data set

Below, the capabilities of `chronospace` are showcased using the
[echinoid data
set](https://github.com/mongiardino/chronospaces_eLife/tree/main/example_files)
described in [Mongiardino Koch et
al. 2022](https://elifesciences.org/articles/72460), comprised of six
posterior distributions of time-calibrated trees obtained using
PhyloBayes (Lartillot et al. 2013). These were all run using the same
constrained topology (as required by the `chronospace` suite of
functions), but varying two methodological decisions:

1)  **Gene subsampling strategy**: three different sets of 100 genes,
    subsampled from a larger phylogenomic data set based on their level
    of clock-likeness or phylogenetic signal, or otherwise at random,
    were assessed.

2)  **Model of molecular evolution**: each gene subsample was also run
    under two models of molecular evolution, the site-homogeneous GTR+G
    and the site-heterogeneous CAT+GTR+G (Lartillot & Philippe 2004).

The code below will download the files into a temporal folder; this
setting will allow us to evaluate whether inferred ages are more
sensitive to the choice of loci or to the choice of model of molecular
evolution.

## Importing and formatting node ages data

The `extract_ages` function imports the input data from a series of
external `.tre` or `.datedist` files (these can be either different runs
of the same analysis or separate analyses), containing trees in Newick
format and stored together in an exclusive folder.

`extract_ages` has two main arguments: `type` - a list of vectors, one
for each factor being tested (in this case two, loci choice and model of
evolution), specifying the group (i.e., factor levels) to which the
chronograms from each file will be assigned to (in this case, three
groups for loci choice and two for model of evolution); and `sample` -
the fixed number of topologies to retain from each file (here we will
sample 500 trees from each file for a total of 3000 trees).

**Note: `extract_ages` will import the files in the order they appear in
the folder; the list provided in `type` must follow this same order**

``` r
#load chronospace
library(chronospace)

#Important! Check files names, compare against type order below
list.files(temp)
#> [1] "clockCATGTR_ln_sample.datedist"  "clockGTR_ln_sample.datedist"    
#> [3] "randomCATGTR_ln_sample.datedist" "randomGTR_ln_sample.datedist"   
#> [5] "signalCATGTR_ln_sample.datedist" "signalGTR_ln_sample.datedist"

#Set type of runs and number of chronograms to be retained. Name the listed 
#elements in type to change factor names:
type <- list(loci = c('clock', 'clock', 'random', 'random', 'signal', 'signal'),
             model = c('CATGTR', 'GTR', 'CATGTR', 'GTR', 'CATGTR', 'GTR'))
sample <- 500

#Import data to R (this might take a minute)
data <- extract_ages(path = temp, type = type, sample = sample)
#> Check that labels are assigned correctly to the input files.
#> If there is an error, modify the order of factors in 'type'
#> or the name of input files, for the two to match.
#> 
#> file = clockCATGTR_ln_sample.datedist | type = clock - CATGTR 
#> file = clockGTR_ln_sample.datedist | type = clock - GTR 
#> file = randomCATGTR_ln_sample.datedist | type = random - CATGTR 
#> file = randomGTR_ln_sample.datedist | type = random - GTR 
#> file = signalCATGTR_ln_sample.datedist | type = signal - CATGTR 
#> file = signalGTR_ln_sample.datedist | type = signal - GTR
```

The function will organize these files into an object of class
`"nodeAges"`, containing 1) a table where each cell represents the
estimated age of a particular node in a particular tree from the
posterior distribution, 2) a data frame where each methodological choice
is coded as a separate factor, and 3) the fixed topology.

``` r
#Print nodeAges object
data
#> Data from 3000 trees with 65 internal nodes (see $ages and $topology),
#> Obtained using 6 methodological pathways (see $factors).
```

## Summarizing chronospaces

The node ages data stored in `"nodeAges"` objects can be summarized
using the `chronospace` function, which uses bgPCA to find the major
directions of variation between groups - i.e., the different
alternatives (levels; e.g., GTR+G vs CAT+GTR+G) associated to a
particular methodological choice (factor; e.g., model of evolution). The
function reports the percentage of variation explained by each factor,
both in relation to the total amount of variation in node ages and for
each individual bgPC extracted. This information can be used to gauge
the relative importance of each choice over node age estimation.

``` r
#Summarize chronospace
cspace <- chronospace(data)
#>                 loci (%) model (%) Unaccounted (%)
#> Total variation  14.0042   5.19656        80.79924
#> _________________________________________________________________________________
#> --- Results for loci (clock/random/signal) ---
#>               loci (%) model (%) Unaccounted (%)
#> bgPC1(17.13%) 10.58553   0.25417         6.28537
#> bgPC2(6.86%)   3.41867   0.09813         3.34235
#> ---------------------------------------------------------------------------------
#> --- Results for model (CATGTR/GTR) ---
#>               loci (%) model (%) Unaccounted (%)
#> bgPC1(11.91%)   0.5823   5.19656         6.12717
#> ---------------------------------------------------------------------------------
#>  * All percentages are relative to the total amount of variation in node ages
```

Chronospaces can be also assessed using visualization tools. The `plot`
function can be used to display scatter plots (when there are three or
more groups) or histograms (when there are only two groups), depicting
the distribution of trees in the space configured by bgPC axes:

``` r
#Plot histogram for model
plot(cspace, factor = 1, output = "ordination")
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

``` r

#Plot scatterplot for loci
plot(cspace, factor = 2, output = "ordination")
```

<img src="man/figures/README-unnamed-chunk-7-2.png" width="100%" />

Another graphical output of using `plot` on a `"nodeAges"` object are
chronograms showing the change in node ages at the extremes of bgPC
axes. Change is represented by color-coding branches according to their
degree of contraction/expansion along the different bgPCA axis to get a
sense of how these different node ages are realized:

``` r
#Plot between-models PC1 extremes 
plot(cspace, factor = 1, output = "extremes")
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-8-2.png" width="100%" />

## Dissecting the effect of methodological choices

`"nodeAges"` objects can also be fed to a number of functions that help
characterizing the impact of methodological choices over particular
aspects of evolutionary history.

`sensitive_nodes` plots the age distribution of the nodes whose ages
vary the most when a particular methodological choice is faced. The user
can either ask for a given number of most affected nodes, or specify a
threshold for the amount of change above which nodes are considered to
be sensitive (arguments `amount_of_change` and `num_clades`,
respectively).

``` r
#Plot the ages distribution of the three nodes most affected by the model of 
#evolution chosen 
sensitive_nodes(echinoid_dates, num_clades = 3, factor = 2)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

`specified_nodes` produces a similar output, but instead of showing the
most affected nodes it displays the ages distribution of a particular
node of interest specified by the user through the argument `tips` (the
function will find the most recent common ancestor between two chosen
tips).

``` r
#Plot the ages distribution of the MRCA of Brissus obesus and Abatus cordatus
#under different loci subsampling strategies
MRCA_Bo_Ac <- specified_node(echinoid_dates, factor = 1,
                             tips = c('Brissus_obesus', 'Abatus_cordatus'))
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

Finally, `ltt_sensitivity` allows exploring the impact of methodological
choices on diversification rates by displaying separate
lineage-through-time curves for each group from a particular factor.
Averaging node ages across subsamples of chronograms can be computed
using either the mean or the median (argument `summary`)

``` r
#Plot lineage through time plots for each model of evolution
ltt_sensitivity(echinoid_dates, factor = 1, summary = "median")
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

## References

Lartillot N., & Philippe H. (2004). *A Bayesian mixture model for
across-site heterogeneities in the amino-acid replacement process*. Mol.
Biol. Evol. 21, 1095–1109. <https://doi.org/10.1093/molbev/msh112>.

Lartillot N., Rodrigue N., Stubbs D., & Richer J. (2013). *PhyloBayes
MPI. Phylogenetic reconstruction with infinite mixtures of profiles in a
parallel environment*. Syst. Biol. 62, 611–615.
<https://doi.org/10.1093/sysbio/syt022>.

Mongiardino Koch N., Thompson J.R., Hiley A.S., McCowin M.F., Armstrong
A.F., Coppard S.E., Aguilera F., Bronstein O., Kroh A., Mooi, R. & Rouse
G.W. (2022). *Phylogenomic analyses of echinoid diversification prompt a
re-evaluation of their fossil record*. Elife, 11, e72460.
<https://doi.org/10.7554/eLife.72460>.
