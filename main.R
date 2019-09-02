library(igraph)
library(ggplot2)

source("./net_metrics.R")

#######   Section 3: Network Population Data   ######
set.seed(123)
# Barabasi-Albert model
ba.pop <- lapply(1:100, function(x) barabasi.game(100, m=5, directed = F))
# Forest Fire model
ff.pop <- lapply(1:100, function(x) forest.fire.game(100, .38, directed = F))
# Stochastic Block model
pm <- matrix(0, nrow=3, ncol=3)
pm[lower.tri(pm)] <- c(p12, p13, p23)
pm <- pm+t(pm)
diag(pm) <- c(p1, p2, p3)
sbm.pop <- lapply(1:100, function(x) sample_sbm(100, pref.matrix = ma))

######   Section 4: Pairwise Evaluation   ######

######   Section 5: Evaluating Variability in Network Populations   ######
