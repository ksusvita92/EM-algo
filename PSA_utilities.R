#
#
# @import dplyr; castor; igraph; network::as.color, mixtools::gammamixEM


require(dplyr)
require(castor)
require(igraph)
require(network)
require(mixtools)
require(reshape2)
require(ggplot2)
require(gridExtra)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Function to get time-phylo
## Input: big tree, cluster (list), last sampled date (to make time-phylo)
## Output: collection of trees (make it as multiphylo), record of infection time
## attr: class "time-phylo"
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
getPhyloPerCluster <- function(tr, clusters, lastSampled) {
  #---- check input completions
  if(class(tr)!="phylo") stop("Tree must be in 'phylo' class")
  if(class(clusters)!="list") stop("clusters must be in 'list' class")
  if(!(all(tr$tip.label %in% unlist(clusters)))) stop("There's at least a tip with no cluster")
  if(!(is.numeric(lastSampled))) stop("Set lastSampled as numeric.")
  
  #---- save clusters and divide into little trees
  df.clust <- NULL
  multi.tree <- tr
  for(cl in names(clusters)) {
    df <- dplyr::data_frame(case=clusters[[cl]],
                     cluster=rep(cl, length(clusters[[cl]])))
    df.clust <- suppressWarnings(bind_rows(df.clust,df))
    little.tr <- get_subtree_with_tips(tr, clusters[[cl]])$subtree
    multi.tree <- c(multi.tree,little.tr)
  }
  df.clust$cluster <- factor(df.clust$cluster, levels = names(clusters))
  names(multi.tree) <- c("all", names(clusters))
  rm(cl,df, little.tr)
  
  #---- get data frame of sampled time
  n<-length(tr$tip.label)
  ed<-tr$edge
  le<-tr$edge.length
  tra<-c(1:n,(2*n-1):(n+1))
  ptree<-matrix(0,2*n-1,3)
  if (n==1) {ptree[1,1]=lastSampled;return(ptree)}
  for (i in 1:nrow(ed)) {
    father<-tra[ed[i,1]]
    son<-tra[ed[i,2]]
    if (ptree[father,2]==0) {ptree[father,2]=son} else {ptree[father,3]=son}
    ptree[son,1]<-le[i]
  }
  todo<-2*n-1
  while (length(todo)>0) {
    t1=todo[1]
    if (ptree[t1,2]==0) {todo=todo[-1];next}
    ptree[ptree[t1,2],1]<-ptree[ptree[t1,2],1]+ptree[t1,1]
    ptree[ptree[t1,3],1]<-ptree[ptree[t1,3],1]+ptree[t1,1]
    todo=c(todo[-1],ptree[t1,2],ptree[t1,3])
  }
  ptree[,1]=ptree[,1]-max(ptree[,1])+lastSampled
  ptree <- suppressWarnings(dplyr::data_frame(case=tr$tip.label,
                             time=ptree[1:n,1]))
  ptree <- ptree %>%
    left_join(df.clust, by="case")
  out <- list(ptree=ptree,
              multiTree=multi.tree)
  class(out) <- "time-phylo"
  return(out)
}
getPhyloFromMetadata <- function(clusters, metadata, tr=NULL) {
  #---- check input completions
  if(!(is.null(tr))) {
    if(class(tr)!="phylo") stop("Tree must be in 'phylo' class")
    if(!(all(tr$tip.label %in% unlist(clusters)))) stop("There's at least a tip with no cluster.")
  }
  if(class(clusters)!="list") stop("clusters must be in 'list' class")
  names(metadata) <- c("case", "time")
  if(class(metadata$case)=="factor") {metadata$name <- with(metadata, as.character(name))}
  if(class(metadata$case)!="character") stop("The cases must in class char or factor.")
  if(class(metadata$time)!="Date") stop("Time must be in Date class.")
  
  #---- get cluster without isolates
  iso.cl_ <- sapply(clusters, function(x) length(x))
  iso.cl_ <- names(which(iso.cl_>1))
  
  #---- save clusters and divide into little trees
  df.clust <- NULL
  multi.tree <- tr
  for(cl in names(clusters)) {
    df <- dplyr::data_frame(case=clusters[[cl]],
                            cluster=rep(cl, length(clusters[[cl]])))
    df.clust <- suppressWarnings(bind_rows(df.clust,df))
    if(!(is.null(multi.tree))) {
      if(length(clusters[[cl]])>1) {
        little.tr <- get_subtree_with_tips(tr, clusters[[cl]])$subtree
        multi.tree <- c(multi.tree,little.tr)
      }
    }
  }
  df.clust$cluster <- factor(df.clust$cluster, levels = names(clusters))
  if(!(is.null(multi.tree))) {names(multi.tree) <- c("all", iso.cl_); rm(little.tr)}
  rm(cl,df)
  
  #---- return as ptree and multiphylo
  ptree <- metadata %>%
    left_join(df.clust, by="case")
  out <- list(ptree=ptree,
              multiTree=multi.tree)
  class(out) <- "time-phylo"
  return(out)
}

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Function to get sampled transmissions
## Input: time-phylo class, # of sampled
## Output: igraph per sampled, transmission net, get sampled serial interval
## attr: class "sampledTT"
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sampleTransTree <- function(tm.phylo, count, replace=F, draw=T) {
  if(class(tm.phylo)!="time-phylo") stop("tm.phylo must be in 'time-phylo' class")
  
  #---- get # of possible transmission trees
  df <- tm.phylo$ptree
  x <- with(df,
            table(cluster))
  n.trans <- prod(factorial(x-1))
  if(count>n.trans) stop(cat("There are ", n.trans, " possible transmission tree to draw."))
  rm(x)
  #
  #---- create possible who infects whom
  dt <- NULL
  clust.nm <- with(df, levels(cluster))
  for(nm in clust.nm) {
    rawdt <- df %>% filter(cluster==nm)
    n <- nrow(rawdt)
    inftor <- rep(rawdt$case, each=n)
    inftee <- rep(rawdt$case, times=n)
    time.i <- rep(rawdt$time, each=n)
    time.j <- rep(rawdt$time, times=n)
    rawdt <- data.frame(i=inftor, 
                        j=inftee,
                        time.i=time.i,
                        time.j=time.j,
                        time.ij=time.j-time.i,
                        cluster=rep(nm,n),
                        stringsAsFactors=F)
    rawdt <- rawdt %>% filter(j!=i & time.ij>=0)
    dt <- bind_rows(dt, rawdt)
  }
  dt$cluster <- factor(dt$cluster, 
                       levels=with(df, levels(cluster)))
  #
  #---- create list to store every possible infectors of a typical case in each cluster
  allTrans_ <- sapply(clust.nm, function(x) NULL) #-- create empty list
  for(nm in clust.nm) {
    dt.clust <- dt %>% filter(cluster==nm)
    x <- unique(dt.clust$j)
    inf.list <- lapply(x, function(x) with(dt.clust, i[j==x])) #-- get all possible infectors of each case
    names(inf.list) <- x
    allTrans_[[nm]] <- inf.list
  }
  #
  #---- get sampled transmission trees from the data
  id_ <- 1 #-- just index
  myGraph <- sapply(paste("draw", 1:count, sep="_"), function(x) NULL)
  while(id_<=count) {
    for(nm in names(allTrans_)) {
      sample_ <- sapply(allTrans_[[nm]], function(x) sample(x,1)) #-- random sample in cluster nm
      df.samp_ <- data.frame(case=sample_,
                             to.case=names(sample_),
                             stringsAsFactors=F)
      myGraph[[id_]] <- bind_rows(myGraph[[id_]], df.samp_)
    }
    if(!replace) { #-- if you want distinct samples
      if(id_>1) {
        for(x in 1:(id_-1)) {
          samedf_ <- identical(myGraph[[id_]], myGraph[[x]]) #-- check if myGraph[[id_]] is identical to previous draws
          if(samedf_) {
            myGraph[[id_]] <- data.frame(case=NULL, to.case=NULL)
            id_ <- id_-1
          }
        }
      }
    } 
    id_ <- id_+1
  }
  rm(id_)
  #
  #---- combine all sampled trans. trees to produce trans. network and eliminate the same edge 
  #     and graph rep
  allEdges <- NULL
  for(id_ in 1:count) {
    allEdges <- bind_rows(allEdges, myGraph[[id_]])
    myGraph[[id_]] <- graph_from_data_frame(myGraph[[id_]], vertices=df)
  }
  allEdges <- allEdges[!(duplicated(allEdges)),]
  allEdges <- allEdges %>% 
    left_join(dt, by=c("case"="i", "to.case"="j")) %>%
    select(case, to.case, time.ij, cluster)
  trans.net <- graph_from_data_frame(allEdges, vertices=df) #-- complete graph
  #
  #---- draw the transmission
  if(draw) {
    plot(trans.net,
         vertex.color=as.color(V(trans.net)$cluster),
         vertex.size=5, 
         vertex.label.cex=1, 
         vertex.label.degree=-pi/2,
         vertex.label.dist=1,
         edge.curve=.3,
         edge.arrow.size=.5,
         layout = layout.fruchterman.reingold(trans.net, niter=10000))
  }
  #
  #---- outputs
  out <- list(count=n.trans, sampled.tt=myGraph, serialInterval=allEdges)
  class(out) <- "sampledTT"
  return(out)
}
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Function to get a data frame of possible serial interval in each cluster
## Input: time-phylo class
## Output: data frame of serial interval
## attr: class "data.frame"
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
getSerialInterval <- function(tm.phylo) {
  if(class(tm.phylo)!="time-phylo") stop("tm.phylo must be in 'time-phylo' class")
  df <- tm.phylo$ptree
  dt <- NULL
  clust.nm <- with(df, levels(cluster))
  for(nm in clust.nm) {
    rawdt <- df %>% filter(cluster==nm)
    n <- nrow(rawdt)
    inftor <- rep(rawdt$case, each=n)
    inftee <- rep(rawdt$case, times=n)
    time.i <- rep(rawdt$time, each=n)
    time.j <- rep(rawdt$time, times=n)
    rawdt <- data.frame(i=inftor, 
                        j=inftee,
                        time.i=time.i,
                        time.j=time.j,
                        time.ij=time.j-time.i,
                        cluster=rep(nm,n),
                        stringsAsFactors=F)
    rawdt <- rawdt %>% filter(j!=i & time.ij>0)
    dt <- bind_rows(dt, rawdt)
  }
  dt$cluster <- factor(dt$cluster, levels=clust.nm)
  return(dt)
}
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Function to get parameter estimation using EM algo and bootstrapping
## Input: time-phylo class, # of trans. sample to draw in each run, bootstrapping iteration, parameter setup
## Output: # of comp., record of boots. params, log-likelihood values in each boots., boots. mean
## attr: class "gammamixBoots"
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
gammamixSetup <- function(lambda0=NULL, alpha0=NULL, beta0=NULL, n.comp=2, tol=1e-8) {
  if(!(is.null(lambda0))) {
    if(class(lambda0)!="numeric") stop("lambda0 must be numeric.")
    if(length(lambda0)!=n.comp) stop("Insufficient # of initial values of lambda0.")
  }
  if(!(is.null(alpha0))) {
    if(class(alpha0)!="numeric") stop("alpha0 must be numeric.")
    if(length(alpha0)!=n.comp) stop("Insufficient # of initial values of alpha0.")
  }
  if(!(is.null(beta0))) {
    if(class(beta0)!="numeric") stop("beta0 must be numeric.")
    if(length(beta0)!=n.comp) stop("Insufficient # of initial values of beta0.")
  }
  #
  #---- outputs
  outputs <- list(n.comp=n.comp,
                  lambda0=lambda0,
                  alpha0=alpha0,
                  beta0=beta0,
                  tol=tol)
  class(outputs) <- "gammamix.control"
  return(outputs)
}
gammamixBoots <- function(tm.phylo, gm.ctrl, ndraw, R=1000) {
  if(class(tm.phylo)!="time-phylo") stop("tm.phylo must be in 'time-phylo' class.")
  if(class(gm.ctrl)!="gammamix.control") stop("gm.ctrl must be in 'gammamix.control' class.")
  #---- how many possible transmission trees?
  ptree <- tm.phylo$ptree
  x <- with(ptree,
            table(cluster))
  n.trans <- prod(factorial(x-1))
  if(ndraw>n.trans) stop(cat("There are ", n.trans, " possible transmission tree to draw."))
  rm(x)
  #
  #---- bootstrapping
  params <- sapply(c("lambda","alpha","beta"), function(x) NULL) 
  r <- 0 
  while(r<=R) {
    si_ <- sampleTransTree(tm.phylo, ndraw, replace=F, draw=F)$serialInterval
    mixEM <- gammamixEM(si_$time.ij, 
                        lambda=gm.ctrl$lambda0,
                        alpha=gm.ctrl$alpha0,
                        beta=gm.ctrl$beta0,
                        k=gm.ctrl$n.comp,
                        epsilon=gm.ctrl$tol,
                        mom.start=F)
    lastiter_ <- length(mixEM$all.loglik)-1
    if(lastiter_<1000) {
      alpha_ <- mixEM$gamma.pars[1,]
      beta_ <- mixEM$gamma.pars[2,]
      lambda_ <- mixEM$lambda
      names(lambda_) <- paste("weight", 1:gm.ctrl$n.comp, sep=".")
      params[["lambda"]] <- bind_rows(params[["lambda"]], lambda_)
      params[["alpha"]] <- bind_rows(params[["alpha"]], alpha_)
      params[["beta"]] <- bind_rows(params[["beta"]], beta_)
      r <- r+1
    }
  }
  #
  #---- exclude outliers in the means
  colMeansq <- function(mat) {
    for(col in 1:ncol(mat)) {
      mat <- as.matrix(mat)
      iqr <- IQR(mat[,col])
      q1 <- quantile(mat[,col], .25)
      q3 <- quantile(mat[,col], .75)
      mild_low <- q1-(1.5*iqr)
      mild_high <- q3+(1.5*iqr)
      mat[which(mat[,col]<=mild_low | mat[,col]>=mild_high)] <- NA
    }
    return(colMeans(mat, na.rm=T))
  }
  mean.params <- lapply(params, function(x) colMeansq(x))
  #
  #---- outputs
  out <- list(n.comp=gm.ctrl$n.comp,
              records=params,
              mean.params=mean.params)
  class(out) <- "gmmixBoots"
  return(out)
}
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Function to plot the bootstrap values of the parameters
## Input: gmmixEM 
## Output: plots
## attr: class "??"
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++








#-----------------------------------------------------
# Define a function to extract clusters
## input: dist. matrix; no. of clusters
## output: cluster size; cluster members
#-----------------------------------------------------
extractCluster <- function(dist, n.clust) {
  if(class(dist)=="dist") { #check the class' variable
    mat <- dist
  } else if(class(dist)=="matrix") {
    if(isSymmetric(dist)) {
      mat <- dist
    } else stop("Message: the matrix is not symmetric.")
  } else stop("Message: incorrect distance matrix.")
  
  getClust <- cutree(hclust(mat), k=n.clust)
  members <- list()
  for(cl in 1:n.clust) members[[cl]] <- names(which(getClust==cl))
  names(members) <- paste("cl", 1:n.clust, sep=".")
  
  return(list(clust.size=sapply(members, FUN = length),
              memberships=members))
}
#

