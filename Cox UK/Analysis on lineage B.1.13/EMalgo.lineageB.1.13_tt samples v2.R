require(dplyr)
#----- Summary and goal --------------------------------------------
#--- //Summary and goal:
#---   - sample X transmission trees for lineage B.1.13 with threshold 0 dSNP
#---   - perform EM-algo in each tree to get a distribution of estimated
#---     shape and scale
#---   - repeat EM-algo with 1 comp., 2 comp., ..., 4 comp. Present this on 
#---     boxplot.
#--- //Steps: * @import EMalgo.lineageB.1.13_dataprep v2.R
#---          * @import PSA_utilities.R, PSA_EM algo.R
#--------------------------------------------------------------- //

#---- dataprep -----------------------------------------------------
#--- //dataprep: * @import EMalgo.lineageB.1.13_dataprep v2.R
#                * @import PSA_utilities.R, PSA_EM algo.R
#----------------------------------------------------------------- //
rootProject <- "/Users/Vita/sfuvault/RProjects/Pairwise Survival Analysis"
setwd(rootProject)
source("PSA_utilities.R")
source("PSA_EM algo.R")

setwd(paste(rootProject, "[Pairwise Survival] COG UK analysis", sep="/"))
source("EMalgo.lineageB.1.13_dataprep v2.R")
#----------------------------------------------------------------- end

#---- parameters ---------------------------------------------------
nsamp_ <- 100 #--- //how many sample to draw?
#----------------------------------------------------------------- //

#---- get sample tt----------------------------------------------------
#--- //goal: * get 100 transmission trees with no repetition
#            * plot some tt in igraph based
#----------------------------------------------------------------- //
data_ <- metaB.1.13 %>% select(sequence_name, sample_date)
tmPhylo_ <- getPhyloFromMetadata(clustThrsh, data_)
tmPhylo1_ <- getPhyloFromMetadata(clustThrsh1, data_)
tmPhylo2_ <- getPhyloFromMetadata(clustThrsh2, data_)
tmPhylo3_ <- getPhyloFromMetadata(clustThrsh, data_)

#set.seed(1234)
sampleTT_ <- sampleTransTree(tmPhylo_, count=nsamp_, draw=F)
sampleTT1_ <- sampleTransTree(tmPhylo1_, count=nsamp_, draw=F)
sampleTT2_ <- sampleTransTree(tmPhylo2_, count=nsamp_, draw=F)
sampleTT3_ <- sampleTransTree(tmPhylo3_, count=nsamp_, draw=F)

#---- apply EM-algo --------------------------------------------------
#--- //goal: * implement EM-algo for each sampled transmission tree
#            * Try different components
#            * store into dataframe
#----------------------------------------------------------------- //
result_ <- NULL
ncomp_ <- 4
si_ <- sampleTT_$serialInterval
si1_ <- sampleTT1_$serialInterval
si2_ <- sampleTT2_$serialInterval
si3_ <- sampleTT3_$serialInterval

temp_ <- data.frame(tt_name=names(sampleTT_$sampled.tt),
                    w1=numeric(nsamp_), w2=numeric(nsamp_),
                    w3=numeric(nsamp_), w4=numeric(nsamp_),
                    alpha=numeric(nsamp_), beta=numeric(nsamp_),
                    error=numeric(nsamp_), niter=numeric(nsamp_),
                    ncomp=numeric(nsamp_), status=character(nsamp_), 
                    stringsAsFactors=F)

for(d_ in c(T,F)) {
  for(comp_ in 1:ncomp_) {
    record_ <- temp_; record_$cutoff <- rep("cutoff: 0 dSNP", nsamp_)
    record1_ <- temp_; record1_$cutoff <- rep("cutoff: 1 dSNP", nsamp_)
    record2_ <- temp_; record2_$cutoff <- rep("cutoff: 2 dSNP", nsamp_)
    record3_ <- temp_; record3_$cutoff <- rep("cutoff: 3 dSNP", nsamp_)
    
    for(nm_ in 1:nsamp_) {
      tt_ <- sampleTT_$sampled.tt[[nm_]]; tt1_ <- sampleTT1_$sampled.tt[[nm_]]
      tt2_ <- sampleTT2_$sampled.tt[[nm_]]; tt3_ <- sampleTT3_$sampled.tt[[nm_]]
      edge_ <- get.edgelist(tt_); edge_ <- as.data.frame(edge_, stringsAsFactors=F)
      edge1_ <- get.edgelist(tt1_); edge1_ <- as.data.frame(edge1_, stringsAsFactors=F)
      edge2_ <- get.edgelist(tt2_); edge2_ <- as.data.frame(edge2_, stringsAsFactors=F)
      edge3_ <- get.edgelist(tt3_); edge3_ <- as.data.frame(edge3_, stringsAsFactors=F)
      names(edge_) <- c("from", "to"); names(edge1_) <- c("from", "to")
      names(edge2_) <- c("from", "to"); names(edge3_) <- c("from", "to")
      #--- //left join with SI
      edge_ <- edge_ %>% left_join(si_, by=c("from"="case", "to"="to.case"))
      edge_$time.ij <- with(edge_, as.numeric(time.ij))
      edge1_ <- edge1_ %>% left_join(si1_, by=c("from"="case", "to"="to.case"))
      edge1_$time.ij <- with(edge1_, as.numeric(time.ij))
      edge2_ <- edge2_ %>% left_join(si2_, by=c("from"="case", "to"="to.case"))
      edge2_$time.ij <- with(edge2_, as.numeric(time.ij))
      edge3_ <- edge3_ %>% left_join(si3_, by=c("from"="case", "to"="to.case"))
      edge3_$time.ij <- with(edge3_, as.numeric(time.ij))
      
      if(d_) {
        est_ <- EMgamma(edge_$time.ij, comp_); est1_ <- EMgamma(edge1_$time.ij, comp_)
        est2_ <- EMgamma(edge2_$time.ij, comp_); est3_ <- EMgamma(edge3_$time.ij, comp_)
        record_$status[nm_] <- "with discrete."; record1_$status[nm_] <- "with discrete."
        record2_$status[nm_] <- "with discrete."; record3_$status[nm_] <- "with discrete."
      } else {
        est_ <- EMgammaND(edge_$time.ij, comp_); est1_ <- EMgammaND(edge1_$time.ij, comp_)
        est2_ <- EMgammaND(edge2_$time.ij, comp_); est3_ <- EMgammaND(edge3_$time.ij, comp_)
        record_$status[nm_] <- "without discrete."; record1_$status[nm_] <- "without discrete."
        record2_$status[nm_] <- "without discrete."; record3_$status[nm_] <- "without discrete."
      }
      
      record_$w1[nm_] <- est_$weights[1]; record_$w2[nm_] <- est_$weights[2]
      record_$w3[nm_] <- est_$weights[3]; record_$w4[nm_] <- est_$weights[4]
      record_$alpha[nm_] <- est_$parameter[1]; record_$beta[nm_] <- est_$parameter[2]
      record_$error[nm_] <- est_$diff; record_$niter[nm_] <- est_$niter 
      record_$ncomp[nm_] <- comp_
      
      record1_$w1[nm_] <- est1_$weights[1]; record1_$w2[nm_] <- est1_$weights[2]
      record1_$w3[nm_] <- est1_$weights[3]; record1_$w4[nm_] <- est1_$weights[4]
      record1_$alpha[nm_] <- est1_$parameter[1]; record1_$beta[nm_] <- est1_$parameter[2]
      record1_$error[nm_] <- est1_$diff; record1_$niter[nm_] <- est1_$niter 
      record1_$ncomp[nm_] <- comp_
      
      record2_$w1[nm_] <- est2_$weights[1]; record2_$w2[nm_] <- est2_$weights[2]
      record2_$w3[nm_] <- est2_$weights[3]; record2_$w4[nm_] <- est2_$weights[4]
      record2_$alpha[nm_] <- est2_$parameter[1]; record2_$beta[nm_] <- est2_$parameter[2]
      record2_$error[nm_] <- est2_$diff; record2_$niter[nm_] <- est2_$niter 
      record2_$ncomp[nm_] <- comp_
      
      record3_$w1[nm_] <- est3_$weights[1]; record3_$w2[nm_] <- est3_$weights[2]
      record3_$w3[nm_] <- est3_$weights[3]; record3_$w4[nm_] <- est3_$weights[4]
      record3_$alpha[nm_] <- est3_$parameter[1]; record3_$beta[nm_] <- est3_$parameter[2]
      record3_$error[nm_] <- est3_$diff; record3_$niter[nm_] <- est3_$niter 
      record3_$ncomp[nm_] <- comp_
    }
    result_ <- bind_rows(result_, record_, record1_, record2_, record3_)
  }
}


#---- output --------------------------------------------------
#--- //goal: * @export [sampleTT] sampleTT_
#            * @export [data.frame] result_
#----------------------------------------------------------------- //
sampledTT <- sampleTT_
sampledTT1 <- sampleTT1_
sampledTT2 <- sampleTT2_
sampledTT3 <- sampleTT3_
EMest <- result_

#--- // clear unnecessary variables
rm(list=grep("_", ls(), value=T))
zz_ <- grep("clust", ls(), value=T)
#--- //remove unnecessary functions
xx_ <- unlist(eapply(.GlobalEnv,class))
yy_ <- names(which(xx_=="function"))
rm(list=c(yy_, zz_))
rm(list=grep("_", ls(), value=T))
rm(metaB.1.13)

cat("@import [data.frame]'EMest': records of EM-algo estimations", "\n")
cat("@import [sampledTT]'sampledTT': 100 sampled transmission trees and serial interval with cutoff=0 dSNP","\n")
cat("@import [sampledTT]'sampledTT1': 100 sampled transmission trees and serial interval with cutoff=1 dSNP","\n")
cat("@import [sampledTT]'sampledTT2': 100 sampled transmission trees and serial interval with cutoff=2 dSNP","\n")
cat("@import [sampledTT]'sampledTT3': 100 sampled transmission trees and serial interval with cutoff=3 dSNP","\n")
