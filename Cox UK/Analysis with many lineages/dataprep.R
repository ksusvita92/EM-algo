#--- //Objective: - filter lineage B, B.1.1, B.1.11, B.1.13, B.1.5, B.1.7, B.2, and B.3
#                 - Choose meaningful cutoffs:
#                     * B: cutoff=0,1,3,4
#                     * B.1.1: cutoff=0,1,2,5
#                     * B.1.11: cutoff=0,4
#                     * B.1.13: cutoff=0,1,3
#                     * B.1.5: cutoff=0,1,3,4
#                     * B.1.7: cutoff=0,1,2
#                     * B.2: cutoff=0,1,2,5
#                     * B.3: cutoff=0,1,3,4
#    //@import: Cox UK data, PSA_clustering.R
#--- //end.

projectRoot <- getwd()
require(dplyr); require(ape)

#---- //@import Cox UK data ----
#     * filter rows that have sequence
#     * filter all wanted lineages
#----------------------------------------------------------------------//
dataRoot_ <- "/Users/Vita/sfuvault/RProjects/Datasets/[Datasets] COG UK"
setwd(dataRoot_)
coxukMeta_ <- read.csv("coxuk_metadata.csv", stringsAsFactors=F)
coxukSeq_ <- read.FASTA("coxuk_sequence.fasta")

coxukMeta_ <- coxukMeta_ %>% filter(sequence_name %in% names(coxukSeq_)) %>%
  filter(lineage %in% c("B","B.1.1","B.1.11","B.1.13","B.1.5","B.1.7","B.2","B.3"))
coxukMeta_$sample_date <- with(coxukMeta_, as.Date(sample_date, "%Y-%m-%d"))
coxukSeq_ <- coxukSeq_[coxukMeta_$sequence_name]
#----------------------------------------------------------------------//

#---- //store each lineage with different cutoffs ----
#     * store each lineage in a list
#     * @import: PSA_clustering.R
#----------------------------------------------------------------------//
lineB <- sapply(paste("cutoff",c(0,1,3,4),sep="="), function(x) NULL)
lineB.1.1 <- sapply(paste("cutoff",c(0,1,2,5),sep="="), function(x) NULL)
lineB.1.11 <- sapply(paste("cutoff",c(0,4),sep="="), function(x) NULL)
lineB.1.13 <- sapply(paste("cutoff",c(0,1,3),sep="="), function(x) NULL)
lineB.1.5 <- sapply(paste("cutoff",c(0,1,3,4),sep="="), function(x) NULL)
lineB.1.7 <- sapply(paste("cutoff",c(0,1,2),sep="="), function(x) NULL)
lineB.2 <- sapply(paste("cutoff",c(0,1,2,5),sep="="), function(x) NULL)
lineB.3 <- sapply(paste("cutoff",c(0,1,3,4),sep="="), function(x) NULL)

clustRoot_ <- "/Users/Vita/sfuvault/RProjects/Pairwise Survival Analysis"
setwd(clustRoot_)
source("PSA_clustering.R")
setwd(projectRoot)  

#--- //get the distance matrix in each lineage
metadtB_ <- coxukMeta_ %>% filter(lineage=="B")
metadtB.1.1_ <- coxukMeta_ %>% filter(lineage=="B.1.1")
metadtB.1.11_ <- coxukMeta_ %>% filter(lineage=="B.1.11")
metadtB.1.13_ <- coxukMeta_ %>% filter(lineage=="B.1.13")
metadtB.1.5_ <- coxukMeta_ %>% filter(lineage=="B.1.5")
metadtB.1.7_ <- coxukMeta_ %>% filter(lineage=="B.1.7")
metadtB.2_ <- coxukMeta_ %>% filter(lineage=="B.2")
metadtB.3_ <- coxukMeta_ %>% filter(lineage=="B.3")

distB_ <- coxukSeq_[metadtB_$sequence_name]; distB_ <- dist.dna(distB_, model="N")
distB.1.1_ <- coxukSeq_[metadtB.1.1_$sequence_name]; distB.1.1_ <- dist.dna(distB.1.1_, model="N")
distB.1.11_ <- coxukSeq_[metadtB.1.11_$sequence_name]; distB.1.11_ <- dist.dna(distB.1.11_, model="N")
distB.1.13_ <- coxukSeq_[metadtB.1.13_$sequence_name]; distB.1.13_ <- dist.dna(distB.1.13_, model="N")
distB.1.5_ <- coxukSeq_[metadtB.1.5_$sequence_name]; distB.1.5_ <- dist.dna(distB.1.5_, model="N")
distB.1.7_ <- coxukSeq_[metadtB.1.7_$sequence_name]; distB.1.7_ <- dist.dna(distB.1.7_, model="N")
distB.2_ <- coxukSeq_[metadtB.2_$sequence_name]; distB.2_ <- dist.dna(distB.2_, model="N")
distB.3_ <- coxukSeq_[metadtB.3_$sequence_name]; distB.3_ <- dist.dna(distB.3_, model="N")

#--- //display plot
mat_ <- matrix(1:8, ncol=4)
layout(mat_, widths=rep.int(1, ncol(mat_)), heights=rep.int(1, nrow(mat_)), respect=F)
hist(distB_, xlab="dSNP", main="Dist. of lineage B")
hist(distB.1.1_, xlab="dSNP", main="Dist. of lineage B.1.1")
hist(distB.1.11_, xlab="dSNP", main="Dist. of lineage B.1.11")
hist(distB.1.13_, xlab="dSNP", main="Dist. of lineage B.1.13")
hist(distB.1.5_, xlab="dSNP", main="Dist. of lineage B.1.5")
hist(distB.1.7_, xlab="dSNP", main="Dist. of lineage B.1.7")
hist(distB.2_, xlab="dSNP", main="Dist. of lineage B.2")
hist(distB.3_, xlab="dSNP", main="Dist. of lineage B.3")

cutoffB_ <- c(0,1,3,4)
cutoffB.1.1_ <- c(0,1,2,5)
cutoffB.1.11_ <- c(0,4)
cutoffB.1.13_ <- c(0,1,3)
cutoffB.1.5_ <- c(0,1,3,4)
cutoffB.1.7_ <- c(0,1,2)
cutoffB.2_ <- c(0,1,2,5)
cutoffB.3_ <- c(0,1,3,4)
  

for(id_ in 1:length(cutoffB_)) {lineB[[id_]] <- getClusters(distB_, cutoffB_[id_])}
for(id_ in 1:length(cutoffB.1.1_)) {lineB.1.1[[id_]] <- getClusters(distB.1.1_, cutoffB.1.1_[id_])}
for(id_ in 1:length(cutoffB.1.11_)) {lineB.1.11[[id_]] <- getClusters(distB.1.11_, cutoffB.1.11_[id_])}
for(id_ in 1:length(cutoffB.1.13_)) {lineB.1.13[[id_]] <- getClusters(distB.1.13_, cutoffB.1.13_[id_])}
for(id_ in 1:length(cutoffB.1.5_)) {lineB.1.5[[id_]] <- getClusters(distB.1.5_, cutoffB.1.5_[id_])}
for(id_ in 1:length(cutoffB.1.7_)) {lineB.1.7[[id_]] <- getClusters(distB.1.7_, cutoffB.1.7_[id_])}
for(id_ in 1:length(cutoffB.2_)) {lineB.2[[id_]] <- getClusters(distB.2_, cutoffB.2_[id_])}
for(id_ in 1:length(cutoffB.3_)) {lineB.3[[id_]] <- getClusters(distB.3_, cutoffB.3_[id_])}

#--- //Show output
cat("lineage B summary:", "\n")
for(id_ in 1:length(cutoffB_)) {
  k_ <- cutoffB_[id_]
  memb_ <- sapply(lineB[[id_]], function(x) length(x))
  nclust_ <- length(which(memb_>1))
  singleton_ <- length(which(memb_==1))
  cat("\t", "cutoff =", cutoffB_[id_], ":", nclust_, "clusters and", singleton_, "singletons.", "\n")
}

cat("lineage B.1.1 summary:", "\n")
for(id_ in 1:length(cutoffB.1.1_)) {
  k_ <- cutoffB.1.1_[id_]
  memb_ <- sapply(lineB.1.1[[id_]], function(x) length(x))
  nclust_ <- length(which(memb_>1))
  singleton_ <- length(which(memb_==1))
  cat("\t", "cutoff =", cutoffB.1.1_[id_], ":", nclust_, "clusters and", singleton_, "singletons.", "\n")
}

cat("lineage B.1.11 summary:", "\n")
for(id_ in 1:length(cutoffB.1.11_)) {
  k_ <- cutoffB.1.11_[id_]
  memb_ <- sapply(lineB.1.11[[id_]], function(x) length(x))
  nclust_ <- length(which(memb_>1))
  singleton_ <- length(which(memb_==1))
  cat("\t", "cutoff =", cutoffB.1.11_[id_], ":", nclust_, "clusters and", singleton_, "singletons.", "\n")
}

cat("lineage B.1.13 summary:", "\n")
for(id_ in 1:length(cutoffB.1.13_)) {
  k_ <- cutoffB.1.13_[id_]
  memb_ <- sapply(lineB.1.13[[id_]], function(x) length(x))
  nclust_ <- length(which(memb_>1))
  singleton_ <- length(which(memb_==1))
  cat("\t", "cutoff =", cutoffB.1.13_[id_], ":", nclust_, "clusters and", singleton_, "singletons.", "\n")
}

cat("lineage B.1.5 summary:", "\n")
for(id_ in 1:length(cutoffB.1.5_)) {
  k_ <- cutoffB.1.5_[id_]
  memb_ <- sapply(lineB.1.5[[id_]], function(x) length(x))
  nclust_ <- length(which(memb_>1))
  singleton_ <- length(which(memb_==1))
  cat("\t", "cutoff =", cutoffB.1.5_[id_], ":", nclust_, "clusters and", singleton_, "singletons.", "\n")
}

cat("lineage B.1.7 summary:", "\n")
for(id_ in 1:length(cutoffB.1.7_)) {
  k_ <- cutoffB.1.7_[id_]
  memb_ <- sapply(lineB.1.7[[id_]], function(x) length(x))
  nclust_ <- length(which(memb_>1))
  singleton_ <- length(which(memb_==1))
  cat("\t", "cutoff =", cutoffB.1.7_[id_], ":", nclust_, "clusters and", singleton_, "singletons.", "\n")
}

cat("lineage B.2 summary:", "\n")
for(id_ in 1:length(cutoffB.2_)) {
  k_ <- cutoffB.2_[id_]
  memb_ <- sapply(lineB.2[[id_]], function(x) length(x))
  nclust_ <- length(which(memb_>1))
  singleton_ <- length(which(memb_==1))
  cat("\t", "cutoff =", cutoffB.2_[id_], ":", nclust_, "clusters and", singleton_, "singletons.", "\n")
}

cat("lineage B.3 summary:", "\n")
for(id_ in 1:length(cutoffB.3_)) {
  k_ <- cutoffB.3_[id_]
  memb_ <- sapply(lineB.3[[id_]], function(x) length(x))
  nclust_ <- length(which(memb_>1))
  singleton_ <- length(which(memb_==1))
  cat("\t", "cutoff =", cutoffB.3_[id_], ":", nclust_, "clusters and", singleton_, "singletons.", "\n")
}
#----------------------------------------------------------------------//

#---- //Outputs ----
coxukMeta <- coxukMeta_
cat("@import clusters from lineages: B, B.1.1, B.1.11, B.1.13, B.1.5, B.1.7, B.2, B.3")
cat("@import Cox UK metadata with above lineages")
rm(list=grep("_", ls(), value=T))
rm(getClusters)
#----------------------------------------------------------------------//
