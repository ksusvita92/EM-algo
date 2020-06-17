#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Functions to produce a data frame of serial interval
## Here we provide a variation of functions to produce time-phylo class.
##    + getPhyloPerCluster requires a tree, cluster-members, and last sampled date
##    + getPhyloFromMetadata requires cluster-members and metadata consisting serial interval
## Input: tree (opt.), cluster-members, last sampled date (opt)
## Output: collection of trees (make it as multiphylo), record of infection time
## attr: class "time-phylo"
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

getPhyloPerCluster <- function(tr, clusters, lastSampled) {
  require(dplyr)
  require(castor) #--- //@import get_subtree_with_tips to subset a tree acc. to the chosen tips
  
  #---- check inputs
  if(class(tr)!="phylo") stop("Tree must be in 'phylo' class")
  if(class(clusters)!="list") stop("clusters must be in 'list' class")
  if(!(all(tr$tip.label %in% unlist(clusters)))) stop("There's at least a tip with no cluster")
  if(!(is.numeric(lastSampled))) stop("Set lastSampled as numeric.")
  
  #---- store clusters and divide into little trees
  df.clust <- NULL
  multi.tree <- tr
  if(is.null(names(clusters))) names(clusters) <- paste("Cl", 1:length(clusters), sep=".")
  for(cl in names(clusters)) {
    df <- data.frame(case=clusters[[cl]],
                     cluster=rep(cl, length(clusters[[cl]])),
                     stingsAsFactors=F)
    df.clust <- suppressWarnings(bind_rows(df.clust,df))
    little.tr <- get_subtree_with_tips(tr, clusters[[cl]])$subtree
    multi.tree <- c(multi.tree,little.tr)
  }
  df.clust$cluster <- factor(df.clust$cluster, levels = names(clusters))
  names(multi.tree) <- c("all", names(clusters))
  rm(cl,df, little.tr)
  
  #---- get data frame of sampled time
  n <- length(tr$tip.label)
  ed <- tr$edge
  le <- tr$edge.length
  tra <- c(1:n,(2*n-1):(n+1))
  ptree <- matrix(0,2*n-1,3)
  if(n==1) {ptree[1,1]=lastSampled;return(ptree)}
  for(i in 1:nrow(ed)) {
    father <- tra[ed[i,1]]
    son <- tra[ed[i,2]]
    if(ptree[father,2]==0) {ptree[father,2] <- son} else {ptree[father,3] <- son}
    ptree[son,1] <- le[i]
  }
  todo <- 2*n-1
  while(length(todo)>0) {
    t1 <- todo[1]
    if(ptree[t1,2]==0) {todo <- todo[-1]; next}
    ptree[ptree[t1,2],1] <- ptree[ptree[t1,2],1]+ptree[t1,1]
    ptree[ptree[t1,3],1] <- ptree[ptree[t1,3],1]+ptree[t1,1]
    todo=c(todo[-1], ptree[t1,2], ptree[t1,3])
  }
  ptree[,1] <- ptree[,1]-max(ptree[,1])+lastSampled
  ptree <- suppressWarnings(data.frame(case=tr$tip.label,
                                       time=ptree[1:n,1]),
                                       stringsAsFactors=F)
  ptree <- ptree %>% left_join(df.clust, by="case")
  out <- list(ptree=ptree,
              multiTree=multi.tree)
  class(out) <- "time-phylo"
  return(out)
}

getPhyloFromMetadata <- function(clusters, metadata, tr=NULL) {
  require(dplyr)
  require(castor) #--- //@import get_subtree_with_tips to subset a tree acc. to the chosen tips

  #---- check input 
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
  iso.cl_ <- sapply(clusters, function(x) length(x)); iso.cl_ <- names(which(iso.cl_>1))
  
  #---- save clusters and divide into little trees
  df.clust <- NULL
  multi.tree <- tr
  if(is.null(names(clusters))) names(clusters) <- paste("Cl", 1:length(clusters), sep=".")
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
  ptree <- metadata %>% left_join(df.clust, by="case")
  out <- list(ptree=ptree,
              multiTree=multi.tree)
  class(out) <- "time-phylo"
  return(out)
}
