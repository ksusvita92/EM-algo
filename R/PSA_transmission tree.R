#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Function to sample transmission trees
## Input: time-phylo class, #of sampled
## Output: igraph per sampled, transmission net (opt), serial interval
## attr: class "sampledTT"
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sampleTransTree <- function(tm.phylo, count, replace=F, draw=T) {
  require(dplyr)
  require(igraph)

  if(class(tm.phylo)!="time-phylo") stop("tm.phylo must be in 'time-phylo' class")
  
  #---- get # of possible transmission trees
  df <- tm.phylo$ptree
  x <- with(df, table(cluster))
  n.trans <- prod(factorial(x-1))
  if(count>n.trans) stop(cat("There are ", n.trans, " possible transmission tree to draw."))
  rm(x)
  
  #---- create possible who-infected-whom and the serial interval
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
  dt$cluster <- factor(dt$cluster, levels=with(df, levels(cluster)))
  
  #---- create a list to store every possible infectors of a typical case in each cluster
  allTrans_ <- sapply(clust.nm, function(x) NULL) 
  for(nm in clust.nm) {
    dt.clust <- dt %>% filter(cluster==nm)
    x <- unique(dt.clust$j)
    inf.list <- lapply(x, function(x) with(dt.clust, i[j==x])) #-- get all possible infectors of each case
    names(inf.list) <- x
    allTrans_[[nm]] <- inf.list
  }
  
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
    if(!replace) { #-- if you want distinct transmission-tree samples
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
