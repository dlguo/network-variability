require(igraph)

# transitivity vector of igraphs list
TransitivitySeries <- function(glist) sapply(glist, transitivity)

# assoratativity of degree of a igraph list
AssortativitySeries <- function(glist) sapply(glist, assortativity_degree)

# local assortativity
LocAssor <- function(gr){
  m1=ecount(gr)
  n1=vcount(gr)
  if(m1==0)
    return(rep(0,n1))
  alpha=numeric(n1)
  tmp=ends(gr,1:m1,names=FALSE)
  d1=degree(gr)
  for(i in 1:m1){
    alpha[tmp[i,1]]=alpha[tmp[i,1]]+((d1[tmp[i,1]]-1)*(d1[tmp[i,2]]-1))/(2*m1)
    alpha[tmp[i,2]]=alpha[tmp[i,2]]+((d1[tmp[i,1]]-1)*(d1[tmp[i,2]]-1))/(2*m1)
  }
  dd=degree_distribution(gr)
  rd=c()
  for(i in 1:(length(dd)-1)){
    rd[i]=(i*dd[(i+1)])/sum((0:(length(dd)-1))*dd)
  }
  rd1=sum((0:(length(rd)-1))*rd)
  beta=(d1*rd1^2)/(2*m1)
  rd2=sum(((0:(length(rd)-1))^2)*rd)
  sig=rd2-rd1^2
  return((alpha-beta)/sig)
}

# local transitivity
LocTrans <- function(gr) transitivity(gr, 'local', isolates = 'zero')

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

DissimSpace <- function(inp.gr, inp.pop){ 	# inp.gr is the input network and inp.pop is the population
  net.count=length(inp.pop)
  inp.dd=degree(inp.gr)
  inp.la=LocAssor(inp.gr)
  inp.lt=LocTrans(inp.gr)
  diss.space=as.data.frame(matrix(0, nrow=net.count, ncol=3))
  colnames(diss.space)=c('Degree', 'Local Assortativity', 'Local Transitivity')
  g1=list()
  for(j in 1:net.count){
    diss.space[j,1]=suppressWarnings(ks.test(inp.dd, degree(inp.pop[[j]])))$statistic
    diss.space[j,2]=suppressWarnings(ks.test(inp.la, LocAssor(inp.pop[[j]])))$statistic
    diss.space[j,3]=suppressWarnings(ks.test(inp.lt, LocTrans(inp.pop[[j]])))$statistic
  }
  return(diss.space)
}

GetDissimModel <- function(inp.pop, simFuncList) {
  size <- length(inp.pop)
  inp.gr <- sample(inp.pop, 1)[[1]]
  dissim.df <- cbind(data.frame(model='real'), as.data.frame(DissimSpace(inp.gr, inp.pop)))
  for (funcName in names(simFuncList)) {
    SimFunc <- simFuncList[[funcName]]
    sim.pop <- SimFunc(inp.gr, size)
    dissim.df <- rbind(dissim.df, cbind(data.frame(model=funcName), as.data.frame(DissimSpace(inp.gr, sim.pop))))
  }
  dissim.df
}