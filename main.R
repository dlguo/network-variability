library(igraph)
library(ergm)
library(ggplot2)
library(GGally)

source("./net_metrics.R")
source("./net_models.R")
#######   Section 3: Network Population Data   ######
set.seed(123)
# Barabasi-Albert model
ba.pop <- lapply(1:100, function(x) barabasi.game(100, m=5, directed = F))
# Forest Fire model
ff.pop <- lapply(1:100, function(x) forest.fire.game(200, .38, directed = F))
# Stochastic Block model
pm <- matrix(0, nrow=3, ncol=3)
pm[lower.tri(pm)] <- c(.07, .07, .07) # The cross-community probability for all communities are 0.07
pm <- pm+t(pm)
diag(pm) <- c(.2, .15, .25) # The inner-community probability for three communities are 0.02, 0.15 and 0.25 respectively.
sbm.pop <- lapply(1:100, function(x) sample_sbm(150, pref.matrix = pm, block.sizes = c(30, 70, 50)))

######   Section 4: Pairwise Evaluation   ######
funcList <- list('Chung-Lu'=SimCL,
                 #'ERGM'=SimERGM,
                 #'ABNG'=SimABNG,
                 'dk2.5'=SimDK)

ba.dissim.df <- GetDissimModel(ba.pop, funcList)
ggpairs(ba.dissim.df, columns = 2:4, mapping = ggplot2::aes(colour=model))
ff.dissim.df <- GetDissimModel(ff.pop, funcList)
ggpairs(ff.dissim.df, columns = 2:4, mapping = ggplot2::aes(colour=model))
sbm.dissim.df <- GetDissimModel(sbm.pop, funcList)
ggpairs(sbm.dissim.df, columns = 2:4, mapping = ggplot2::aes(colour=model))

######   Section 5: Evaluating Variability in Network Populations   ######
