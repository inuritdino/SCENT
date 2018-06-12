### CompSRanaPRL.R
### Author: Andrew Teschendorff
### Date: 12 Jun 2018
### This software is released under a GPL v3 licence.

### AUXILIARY FUNCTIONS ##################################################

### CompS: Computes local entropy of a node with stochastic vector p.v
CompS <- function(p.v){
  tmp.idx <- which(p.v>0);
  S <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )
  return(S);
}

### CompNS: Computes the normalised local entropy of a node with stochastic vector p.v
CompNS <- function(p.v){
  tmp.idx <- which(p.v>0);
  if(length(tmp.idx)>1){
    S <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )/log(length(tmp.idx));
  }
  else { ### one degree nodes have zero entropy, avoid singularity.
    S <- 0;
  }
  return(S);
}

### betaFn: Linear function defining gene expression values of a given sample (exp.v) into low, intermediate and high expression levels according to two thresholds specified by vector alpha.v. Gene expression is assumed to come from Illumina or Affy arrays. Returns expression values normalised between 0 (low expression) and 1 (high expression) with intermediate values taking a linear function between 0 and 1.

betaFn <- function(exp.v,alpha.v){
  high.idx <- which(exp.v>=alpha.v[2]);
  low.idx <-  which(exp.v<=alpha.v[1]);
  mid.idx <- setdiff(1:length(exp.v),c(high.idx,low.idx));
  out.v <- exp.v;
  out.v[high.idx] <- 1;
  out.v[low.idx] <- 0;
  out.v[mid.idx] <- (exp.v[mid.idx]-alpha.v[1])/(alpha.v[2]-alpha.v[1]);
  return(out.v);
}

### CompMaxSR: Computes the maximum entropy rate of a network with adjacency matrix adj.m (assumed connected).

CompMaxSR <- function(adj.m){

require(igraph);
# find right eigenvector of adjacency matrix
fa <- function(x,extra=NULL) {
   as.vector(adj.m %*% x)
}
ap.o <- arpack(fa, options=list(n=nrow(adj.m),nev=1,which="LM"),sym=TRUE);
v <- ap.o$vectors;
lambda <- ap.o$values;
maxSR <- log(lambda); ### maximum entropy
return(maxSR);

}


##############################################################################

### MAIN FUNCTION

### CompSR: This function computes the entropy rate (SR) of a sample specified by a gene expression profile exp.v in a network specified by an adjacency matrix adj.m. The latter is assumed connected.
### INPUT:
### adj.m: an adjacency matrix representing a connected PPI network, with rownames/colnames of the matrix annotated to Entrez gene IDs.
### exp.v: a genome-wide expression vector of a sample with names labeling the genes (also annotated to entrez gene IDs) of same length as rows of adj.m.
### k.v: the degree distribution of adj.m, i.e. k.v=rowSums(adj.m)
### local: logical, if TRUE also compute the local Shannon (normalised) entropies.
### method: scalar specifying method to compute stochastic matrix: 1=stochastic vector of a node is independent of the node's gene expression value, 2=stochastic vector of a node is dependent on the node's gene expression value through use of the betaFn. Default is 1.
### quants.v: optionally, a vector of length 2 specifying the quantiles for defining low and high expression (only needed if method=2 is used).
### maxSR: optionally, the maximum entropy rate of the network.

### OUTPUT: a list of four objects
### sr: the entropy rate of the sample (normalised between 0 and 1 if maxSR was provided).
### inv: the stationary distribution of the sample.
### s: the unnormlised local entropies of the sample.
### ns: the normalised local entropies of the sample.

CompSR <- function(exp.v,adj.m,k.v,local=TRUE,method=1,quants.v=c(0.1,0.9),maxSR=NULL){

  require(igraph);
  ### compute stochastic matrix
  p.m <- matrix(0,nrow=length(exp.v),ncol=length(exp.v));
  rownames(p.m) <- rownames(adj.m);
  colnames(p.m) <- rownames(adj.m);

  if(method==1){
   for(g in 1:nrow(adj.m)){
    nn.idx <- which(adj.m[g,]==1);
    term2.v <- exp.v[nn.idx]/sum(exp.v[nn.idx]);
    p.m[g,nn.idx] <- term2.v;
   }
  }
  else if(method==2){
   alpha.v <- quantile(as.vector(exp.v),probs=quants.v);
   beta.v <- betaFn(exp.v,alpha.v);
   for(g in 1:nrow(adj.m)){
    nn.idx <- which(adj.m[g,]==1);
    term2.v <- exp.v[nn.idx]/sum(exp.v[nn.idx]);
    p.m[g,nn.idx] <- (1-beta.v[g])/k.v[g] + beta.v[g]*term2.v
   }
  }
 # print("Computed stochastic matrix");

  fp <- function(x,extra=NULL) {
   as.vector(p.m %*% x)
  }
  fpt <- function(x,extra=NULL) {
   as.vector(t(p.m) %*% x)
  }
  ap.o <- arpack(fpt, options=list(n=nrow(p.m),nev=1,which="LM"),sym=FALSE);
  invP.v <- abs(as.numeric(ap.o$vectors));
  invP.v <- invP.v/sum(invP.v);
  S.v <- apply(p.m,1,CompS);
  SR <- sum(invP.v*S.v);
  if(is.null(maxSR)==FALSE){## if provided then normalise relative to maxSR
    SR <- SR/maxSR;
  }
  if(local){
   NS.v <- apply(p.m,1,CompNS);
  }
  else {
   NS.v <- NULL;
  }

  return(list(sr=SR,inv=invP.v,s=S.v,ns=NS.v));
}

###################################################################################


### MAIN FUNCTION

### CompSRana: This function computes the entropy rate (SR) of a sample specified by a gene expression profile exp.v in a network specified by an adjacency matrix adj.m. The latter is assumed connected.
### INPUT:
### adj.m: an adjacency matrix representing a connected PPI network, with rownames/colnames of the matrix annotated to Entrez gene IDs.
### exp.v: a genome-wide expression vector of a sample with names labeling the genes (also annotated to entrez gene IDs) of same length as rows of adj.m.
### k.v: the degree distribution of adj.m, i.e. k.v=rowSums(adj.m)
### local: logical, if TRUE also compute the local Shannon (normalised) entropies.
### method: scalar specifying method to compute stochastic matrix: 1=stochastic vector of a node is independent of the node's gene expression value, 2=stochastic vector of a node is dependent on the node's gene expression value through use of the betaFn. Default is 1.
### quants.v: optionally, a vector of length 2 specifying the quantiles for defining low and high expression (only needed if method=2 is used).
### maxSR: optionally, the maximum entropy rate of the network.

### OUTPUT: a list of four objects
### sr: the entropy rate of the sample (normalised between 0 and 1 if maxSR was provided).
### inv: the stationary distribution of the sample.
### s: the unnormlised local entropies of the sample.
### ns: the normalised local entropies of the sample.


CompSRana <- function(exp.v,adj.m,local=TRUE,maxSR=NULL){

  ### compute outgoing flux around each node
  sumexp.v <- as.vector(adj.m %*% matrix(exp.v,ncol=1));
  invP.v <- exp.v*sumexp.v;
  nf <- sum(invP.v);
  invP.v <- invP.v/nf;
  p.m <- t(t(adj.m)*exp.v)/sumexp.v;
  S.v <- apply(p.m,1,CompS);
  SR <- sum(invP.v*S.v);
  if(is.null(maxSR)==FALSE){## if provided then normalise relative to maxSR
   SR <- SR/maxSR;
  }
  if(local){
   NS.v <- apply(p.m,1,CompNS);
  }
  else {
   NS.v <- NULL;
  }
  return(list(sr=SR,inv=invP.v,s=S.v,ns=NS.v));
}


CompSRanaPRL <- function(idx,exp.m,adj.m,local=TRUE,maxSR=NULL){

  ### compute outgoing flux around each node
  exp.v <- exp.m[,idx];
  sumexp.v <- as.vector(adj.m %*% matrix(exp.v,ncol=1));
  invP.v <- exp.v*sumexp.v;
  nf <- sum(invP.v);
  invP.v <- invP.v/nf;
  p.m <- t(t(adj.m)*exp.v)/sumexp.v;
  S.v <- apply(p.m,1,CompS);
  SR <- sum(invP.v*S.v);
  if(is.null(maxSR)==FALSE){## if provided then normalise relative to maxSR
   SR <- SR/maxSR;
  }
  if(local){
   NS.v <- apply(p.m,1,CompNS);
  }
  else {
   NS.v <- NULL;
  }
  return(list(sr=SR,inv=invP.v,s=S.v,ns=NS.v));
}


#### Mixing rate functions
CompMixRate <- function(exp.v,adj.m,nLV=5){
  sumexp.v <- as.vector(adj.m %*% matrix(exp.v,ncol=1));
  invP.v <- exp.v*sumexp.v;
  nf <- sum(invP.v);
  invP.v <- invP.v/nf;
  p.m <- t(t(adj.m)*exp.v)/sumexp.v;
  require(igraph);
  # find right eigenvector of stochastic matrix
  fp <- function(x,extra=NULL) {
   as.vector(p.m %*% x)
  }
  ap.o <- arpack(fp, options=list(n=nrow(p.m),nev=2,ncv=nLV,which="LM"),sym=FALSE);
  slem <- abs(as.numeric(ap.o$value)[2]);
  muR <- log(1/slem);
  return(muR);

}

CompMixRatePRL <- function(idx,exp.m,adj.m,nLV=5){
  exp.v <- exp.m[,idx];
  sumexp.v <- as.vector(adj.m %*% matrix(exp.v,ncol=1));
  invP.v <- exp.v*sumexp.v;
  nf <- sum(invP.v);
  invP.v <- invP.v/nf;
  p.m <- t(t(adj.m)*exp.v)/sumexp.v;
  require(igraph);
  # find right eigenvector of stochastic matrix
  fp <- function(x,extra=NULL) {
   as.vector(p.m %*% x)
  }
  ap.o <- arpack(fp, options=list(n=nrow(p.m),nev=2,ncv=nLV,which="LM",tol=0.0001),sym=FALSE);
  slem <- abs(as.numeric(ap.o$value)[2]);
  muR <- log(1/slem);
  return(muR);

}
