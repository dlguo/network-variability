require(igraph)
require(ergm)

# Chung-Lu model
SimCL <- function(inp.gr, size) lapply(1:size, function(x) static.fitness.game(ecount(inp.gr),degree(inp.gr)))

# Exponential Random Graph Model
SimERGM <- function(inp.gr, size) {
  adj=as.matrix(get.adjacency(inp.gr))
  net=network(adj, matrix.type="adjacency", directed=FALSE)
  fit1 <- ergm(net ~ edges + gwdegree(1,fixed=FALSE)+gwesp(1,fixed=FALSE)+gwdsp(0,fixed=FALSE), control=control.ergm(MCMLE.maxit=50))    #MCMLE.maxit can be set appropriately
  lapply(1:size, function(x) graph.adjacency(as.matrix(simulate(fit1)), 'undirected'))
}

SimDK <- function(inp.gr, size) {
  if(!file.exists('RandNetGen/RandNetGen')) system('./installDK.sh')
  write_graph(inp.gr, './inpgr', 'edgelist')
  DK.pop <- list()
  for (i in 1:size) {
    system('./RandNetGen/RandNetGen -net inpgr -dk 2.5')
    DK.pop[[i]] <- simplify(graph_from_edgelist(as.matrix(read.table('./dk2.5_inpgr'))+1, directed = F))
    file.remove(c('./dk2.5_inpgr', './E_vs_T.dat'))
  }
  file.remove('./inpgr')
  DK.pop
}
