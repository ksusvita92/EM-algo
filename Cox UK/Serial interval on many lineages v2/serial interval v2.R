#---- //Some notes ----
# - The difference between "serial interval.R" is that in here we put restriction on
#   estimated shape, such that alpha>1.
# - We choose the cutoff=1 for all lineages, but lineage B.1.11. 
#   We choose cutoff=0 for lineage B.1.11.
#----------------------------------------------------------------------------------------//

source("dataprep.R")

setwd("/Users/Vita/sfuvault/RProjects/Pairwise Survival Analysis")
source("PSA_utilities.R")
source("PSA_EM algo with restriction.R")
setwd(projectRoot)

#---- //Sample n transmission trees ----
set.seed(1234)
nsamp_ <- 500
ncomp_ <- 1
ttSamples <- sapply(paste("lineage", c("B","B.1.1","B.1.11","B.1.13","B.1.5","B.1.7","B.2","B.3"), sep=""), function(x) NULL)

for(nm_ in names(ttSamples)) {
  if(nm_=="lineageB") {
    metadt_ <- coxukMeta %>% filter(lineage=="B") %>% select(sequence_name, sample_date)
    tmPhylo_ <- getPhyloFromMetadata(lineB[[2]], metadt_)
  } else if(nm_=="lineageB.1.1") {
    metadt_ <- coxukMeta %>% filter(lineage=="B.1.1") %>% select(sequence_name, sample_date)
    tmPhylo_ <- getPhyloFromMetadata(lineB.1.1[[2]], metadt_)
  } else if(nm_=="lineageB.1.11") {
    metadt_ <- coxukMeta %>% filter(lineage=="B.1.11") %>% select(sequence_name, sample_date)
    tmPhylo_ <- getPhyloFromMetadata(lineB.1.11[[2]], metadt_)
  } else if(nm_=="lineageB.1.13") {
    metadt_ <- coxukMeta %>% filter(lineage=="B.1.13") %>% select(sequence_name, sample_date)
    tmPhylo_ <- getPhyloFromMetadata(lineB.1.13[[2]], metadt_)
  } else if(nm_=="lineageB.1.5") {
    metadt_ <- coxukMeta %>% filter(lineage=="B.1.5") %>% select(sequence_name, sample_date)
    tmPhylo_ <- getPhyloFromMetadata(lineB.1.5[[2]], metadt_)
  } else if(nm_=="lineageB.1.7") {
    metadt_ <- coxukMeta %>% filter(lineage=="B.1.7") %>% select(sequence_name, sample_date)
    tmPhylo_ <- getPhyloFromMetadata(lineB.1.7[[2]], metadt_)
  } else if(nm_=="lineageB.2") {
    metadt_ <- coxukMeta %>% filter(lineage=="B.2") %>% select(sequence_name, sample_date)
    tmPhylo_ <- getPhyloFromMetadata(lineB.2[[2]], metadt_)
  } else if(nm_=="lineageB.3") {
    metadt_ <- coxukMeta %>% filter(lineage=="B.3") %>% select(sequence_name, sample_date)
    tmPhylo_ <- getPhyloFromMetadata(lineB.3[[2]], metadt_)
  }
  ttSamples[[nm_]] <- sampleTransTree(tmPhylo_, nsamp_, draw=F)
}
#-------------------------------------------------------------------------//

#---- //Perform EM-algo ----
result <- NULL
for(nm_ in names(ttSamples)) {
  sampleTT_ <- ttSamples[[nm_]]
  si_ <- sampleTT_$serialInterval
  temp_ <- data.frame(ttName=c(names(sampleTT_$sampled.tt),"mean_"),
                      w1=numeric(nsamp_+1), w2=numeric(nsamp_+1),
                      w3=numeric(nsamp_+1), w4=numeric(nsamp_+1),
                      alpha=numeric(nsamp_+1), beta=numeric(nsamp_+1),
                      mean=numeric(nsamp_+1), lineage=rep(nm_, nsamp_+1),
                      status=c(rep("boot",nsamp_),"mean"), stringsAsFactors=F)
  for(id_ in 1:nsamp_) {
    ttree_ <- sampleTT_$sampled.tt[[id_]]
    edge_ <- igraph::get.edgelist(ttree_); edge_ <- as.data.frame(edge_, stringsAsFactors=F)
    names(edge_) <- c("from", "to")
    edge_ <- edge_ %>% left_join(si_, by=c("from"="case", "to"="to.case"))
    edge_$time.ij <- with(edge_, as.numeric(time.ij))
    
    est_ <- EMgammaND(edge_$time.ij, ncomp_)
    if(est_$niter<1000) {
      temp_$w1[id_] <- est_$weights[1]; temp_$w2[id_] <- est_$weights[2]; temp_$w3[id_] <- est_$weights[3]; temp_$w4[id_] <- est_$weights[4]
      temp_$alpha[id_] <- est_$parameter[1]; temp_$beta[id_] <- est_$parameter[2]
      temp_$mean[id_] <- prod(est_$parameter)
    } else {
      temp_$w1[id_] <- NA; temp_$w2[id_] <- NA; temp_$w3[id_] <- NA; temp_$w4[id_] <- NA
      temp_$alpha[id_] <- NA; temp_$beta[id_] <- NA
      temp_$mean[id_] <- NA
    }
    
  }
  temp_$w1[nsamp_+1] <- with(temp_ %>% filter(ttName!="mean_"), mean(w1, na.rm=T))
  temp_$w2[nsamp_+1] <- with(temp_ %>% filter(ttName!="mean_"), mean(w2, na.rm=T))
  temp_$w3[nsamp_+1] <- with(temp_ %>% filter(ttName!="mean_"), mean(w3, na.rm=T))
  temp_$alpha[nsamp_+1] <- with(temp_ %>% filter(ttName!="mean_"), mean(alpha, na.rm=T))
  temp_$beta[nsamp_+1] <- with(temp_ %>% filter(ttName!="mean_"), mean(beta, na.rm=T))
  temp_$mean[nsamp_+1] <- with(temp_[nsamp_+1,], alpha*beta)
  result <- bind_rows(result, temp_)
}
rm(list=grep("_",ls(),value=T))
rm(list=grep("lineB", ls(), value=T))
#-------------------------------------------------------------------------//


#---- //Plot the distribution of serial interval ----
data <- NULL
n_ <- 100
x_ <- seq(.01, 10, length.out=n_)

for(id_ in 1:nrow(result)) {
  alpha_ <- result$alpha[id_]
  beta_ <- result$beta[id_]
  temp_ <- data.frame(ttName=rep(result$ttName[id_], n_),
                      x=x_, fx=dgamma(x_, alpha_, 1/beta_),
                      lineage=rep(result$lineage[id_], n_),
                      status=rep(result$status[id_], n_),
                      stringsAsFactors=F)
  data <- bind_rows(data, temp_)
}
rm(list=grep("_",ls(),value=T))

plot1 <- ggplot(data=(data %>% filter(status=="boot")), aes(x=x, y=fx, group=ttName)) + 
  geom_line(col="gray", lwd=1) + 
  geom_line(data=(data %>% filter(status=="mean")), col="blue") +
  facet_wrap(~lineage, nrow=2) +
  xlab("serial interval") + ylab("density") + 
  labs(title="Distibution of serial interval", subtitle="with ncomp=1 and 500 bootstrap")
plot2 <- ggplot(data %>% filter(status=="mean")) +
  geom_line(aes(x=x, y=fx, col=lineage)) +
  xlab("serial interval") + ylab("density") +
  labs(title="Distibution of serial interval", subtitle="with ncomp=1 and 500 bootstrap")

#---- //Plot the distribution of mean of serial interval from bootstrapping ----
plot3 <- ggplot(result %>% filter(status=="boot")) +
  geom_histogram(aes(x=mean), alpha=.4, fill="red") +
  facet_wrap(~lineage, nrow=2) + theme(legend.position="") +
  labs(title="Dist. of mean serial interval", subtitle="with ncomp=1 and from 500 bootstrapping") +
  xlab("mean serial interval")

result %>% filter(status=="mean") %>% select(lineage, mean)

