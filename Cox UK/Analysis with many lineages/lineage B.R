#--- //Do analysis for lineage B
#    Objective: - sample n transmission trees with no repetition
#               - perform EM-algo with ncomp=1,2,3,4 in each tree
#               - get the estimation of alpha, beta, and mixture mean
#    @import: PSA_utilities.R, PSA_EM algo.R, igraph::get.edgelist
#----------------------------------------------------------------------//

setwd("/Users/Vita/sfuvault/RProjects/Pairwise Survival Analysis")
source("PSA_utilities.R")
source("PSA_EM algo.R")
setwd("/Users/Vita/sfuvault/RProjects/Pairwise Survival Analysis/[Pairwise Survival] COG UK analysis/With many lineages")
source("dataprep.R")
rm(list=grep("B.", ls(), value=T))
metadtB_ <- coxukMeta %>% filter(lineage=="B") %>% select(sequence_name, sample_date)
nsamp_ <- 100 #--- //how many trees to sample?

#---- //sample nsamp_ trans. trees ----
#     * Store into a list
#----------------------------------------------------------------------//
ttSample_ <- sapply(names(lineB), function(x) NULL)
for(nm_ in names(ttSample_)) {
  tmPhylo_<- getPhyloFromMetadata(lineB[[nm_]], metadtB_)
  ttSample_[[nm_]] <- sampleTransTree(tmPhylo_, nsamp_, draw=F)
}
#----------------------------------------------------------------------//


#---- //perform EM-algo in each trans. tree and each cutoff
#     * get dataframe as a summary
#     * we perform EM-algo with ncomp=1,2,3,4
#----------------------------------------------------------------------//
set.seed(12345)
ncomp_ <- 4
result_ <- sapply(names(ttSample_), function(x) NULL)
obsSI_ <- sapply(names(ttSample_), function(x) NULL)

for(nm_ in names(result_)) {
  transSamp_ <- ttSample_[[nm_]]
  si_ <- transSamp_$serialInterval
  for(comp_ in 1:ncomp_) {
    record_ <- data.frame(ttName=names(transSamp_$sampled.tt),
                          ncomp=rep(paste("ncomp", comp_, sep="="), nsamp_),
                          w1=numeric(nsamp_), w2=numeric(nsamp_),
                          w3=numeric(nsamp_), w4=numeric(nsamp_),
                          alpha=numeric(nsamp_), beta=numeric(nsamp_),
                          mixtureMean=numeric(nsamp_),
                          stringsAsFactors=F)
    for(id_ in 1:nsamp_) {
       ttree_ <- transSamp_$sampled.tt[[id_]]
       edge_ <- igraph::get.edgelist(ttree_); edge_ <- as.data.frame(edge_, stringsAsFactors=F)
       names(edge_) <- c("from", "to")
       edge_ <- edge_ %>% left_join(si_, by=c("from"="case", "to"="to.case"))
       edge_$time.ij <- with(edge_, as.numeric(time.ij))
       
       est_ <- EMgammaND(edge_$time.ij, comp_)
       record_$w1[id_] <- est_$weights[1]; record_$w2[id_] <- est_$weights[2]; record_$w3[id_] <- est_$weights[3]; record_$w4[id_] <- est_$weights[4]
       record_$alpha[id_] <- est_$parameter[1]; record_$beta[id_] <- est_$parameter[2]
       record_$mixtureMean[id_] <- sum(est_$weights*(1:4)*prod(est_$parameter))
    }
    result_[[nm_]] <- bind_rows(result_[[nm_]], record_)
  }
  obsSI_[[nm_]] <- si_
}
#----------------------------------------------------------------------//


#---- //Get a tidy output ----
#     * bind all observation with different cutoffs
#----------------------------------------------------------------------//
record_ <- NULL

for(nm_ in names(result_)) {
  temp_ <- result_[[nm_]]
  temp_$cutoff <- rep(nm_, nrow(temp_))
  record_ <- bind_rows(record_, temp_)
}
summary_ <- record_ %>% group_by(ncomp, cutoff) %>%
  summarise(meanSI=mean(mixtureMean), min=min(mixtureMean), max=max(mixtureMean))
#----------------------------------------------------------------------//


#---- //outputs ----
record <- record_
rm(list=grep("_", ls(), value=T))
rm(EMgamma, EMgammaND, extractCluster, gammamixSetup, getPhyloFromMetadata, getPhyloPerCluster, getSerialInterval, sampleTransTree)
#----------------------------------------------------------------------//
