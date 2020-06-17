#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Function to extract clusters from distance matrix and given threshold.
## Input: 'dist' matrix and 'double' threshold
## Output: a list of formed clusters
## attr: class 'list'
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
require(dplyr)

getClusters <- function(dist, thrsh) {
  #--- //Func. to change dist to data frame
  distToDF = function(inDist, val.name) {
    if (class(inDist) != "dist") stop("wrong input type")
    A <- attr(inDist, "Size")
    B <- if (is.null(attr(inDist, "Labels"))) paste("case",sequence(A)) else attr(inDist, "Labels")
    if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
    if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
    
    df = data.frame(
      row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
      col = rep(B[-length(B)], (length(B)-1):1),
      value = as.vector(inDist),
      stringsAsFactors=F)
    names(df) = c("row", "col", val.name)
    return(df)
  }
  
  FulldistDF_ <- distToDF(dist, "distance")
  distDF_ <- FulldistDF_ %>% filter(distance<=thrsh) #--- //subset DF <= threshold
  caseName <- with(distDF_, unique(c(row,col)))
  Cluster_ <- list()
  member_ <- NULL
  case <- caseName[1]; caseIndex <- caseName[1]
  id_ <- 1
  
  while(case %in% caseIndex & !(is.na(case))) {
    pairCol_ <- which(distDF_$row==case); pairRow_ <- which(distDF_$col==case)
    pairCol_ <- distDF_$col[pairCol_]; pairRow_ <- distDF_$row[pairRow_]
    vv_ <- c(case, pairCol_, pairRow_) #--- //get all pairwise cases
    
    if(length(which(vv_ %in% member_))==0) { #--- //check if some vv_ has been clustered
      clust_ <- list(vv_) #--- //define new cluster
      names(clust_) <- paste("Cl", id_, sep=".")
      Cluster_ <- append(Cluster_, clust_)
      id_ <- id_+1
    } else {
      where_ <- sapply(Cluster_, function(x) vv_ %in% x)
      col_ <- 1
      temp_ <- F
      while(col_<=ncol(where_) & temp_==F) {
        check_ <- any(where_[,col_])
        if(check_) {
          temp_ <- T 
          name_ <- colnames(where_)[col_]
        } else {col_ <- col_+1}
      }
      Cluster_[[name_]] <- c(Cluster_[[name_]], vv_)
      Cluster_[[name_]] <- unique(Cluster_[[name_]])
    }
    member_ <- unlist(Cluster_, F, F)
    caseIndex <- setdiff(caseName, member_)
    case <- caseIndex[1]
  }
  
  #--- //now include all isolates in different cluster
  allCases <- with(FulldistDF_, unique(c(row,col)))
  unclust_ <- setdiff(allCases, caseName)
  if(length(unclust_)!=0) {
    for(each_ in 1:length(unclust_)) {
      #id_:(length(unclust_)+id_-1)
      clust_ <- list(unclust_[each_])
      names(clust_) <- paste("Cl", id_, sep=".")
      Cluster_ <- append(Cluster_, clust_)
      id_ <- id_+1
    }
  }
  return(Cluster_)
}
