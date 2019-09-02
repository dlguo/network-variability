require(igraph)

# transitivity vector of igraphs list
TransitivitySeries <- function(glist) sapply(glist, transitivity)

# assoratativity of degree of a igraph list
AssortativitySeries <- function(glist) sapply(glist, assortativity_degree)

# entropy of a probability vector
entropyProb <- function(prob, base=length(prob)) {
  if (length(prob) == 1) return(0)
  else return(-sum(sapply(prob, function(p) if (p==0) 0 else p * log(p, base))))
}

# edge existence entropy
edgeEntropy <- function(adjmat) {
  rList <- list()
  N <- dim(adjmat)[3]
  n <- dim(adjmat)[1]
  rList[['N']] <- N
  probmat <- rowSums(adjmat, dims=2)/N
  inputmat <- cbind(c(probmat), 1-c(probmat))
  entropy_mat <- matrix(apply(inputmat, MARGIN = 1, function(x) entropyProb(x,base=2)), nrow=n)
  rList[['entropy_mat']]<- entropy_mat
  rList[['entropy_mean']] <- mean(entropy_mat)
  rList[['entropy_mean_offdiag']] <- mean(entropy_mat[lower.tri(entropy_mat, diag = F)])
  rList
}

# geodesic entropy
GeodesicEntropy <- function(adjmat, rm.neighbor=F) {
  rList <- list()
  N <- dim(adjmat)[3]
  n <- dim(adjmat)[1]
  rList[['N']] <- N
  distmat <- array(dim=c(n,n,N))
  for (i in 1:N) {
    distmat[,,i] <- distances(graph_from_adjacency_matrix(adjmat[,,i], mode = "undirected", diag = F))
  }
  entropy_mat <- array(0, dim=c(n,n))
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (rm.neighbor) {
        td <- table(distmat[i,j,])
        td <- td[-which(names(td)=='1')]
        dist_freq <- td/sum(td)
      } else {
        dist_freq <- table(distmat[i,j,])/N
      }
      if (length(dist_freq)>1) entropy_mat[i,j] <- entropyProb(dist_freq)
    }
  }
  entropy_mat <- entropy_mat + t(entropy_mat)
  rList[['entropy_mat']]<- entropy_mat
  rList[['entropy_mean']] <- mean(entropy_mat)
  rList[['entropy_mean_offdiag']] <- mean(entropy_mat[lower.tri(entropy_mat, diag = F)])
  rList
}