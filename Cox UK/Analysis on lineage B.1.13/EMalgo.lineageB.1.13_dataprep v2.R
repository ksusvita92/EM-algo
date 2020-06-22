#----- Summary and goal --------------------------------------------
#--- //Summary and goal:
#---   - extract lineage B.1.13 of COG UK
#---   - get its cluster according to the choosen cutoff
#--- //Steps: * @import PSA_clustering.R
#--------------------------------------------------------------- //

#---- parameters ---------------------------------------------------
thrsh_ <- 0 #--- //based on its median
#--------------------------------------------------------------- end


#---- dataprep -----------------------------------------------------
#--- //dataprep: * @import clustering R-script
#                * @import COG UK dataset and @extract lineage B.1.13
#----------------------------------------------------------------- //
rootProject <- "/Users/Vita/sfuvault/RProjects"
setwd(paste(rootProject, "Pairwise Survival Analysis", sep="/"))
source("PSA_clustering.R")

setwd(paste(rootProject,"Datasets/[Datasets] COG UK", sep="/"))
metadata_ <- read.csv("coxuk_metadata.csv", stringsAsFactors=F)
metadata_$sample_date <- with(metadata_, as.Date(sample_date, "%Y-%m-%d"))
require(ape)
sequence_ <- read.FASTA("coxuk_sequence.fasta")

withSeq_ <- metadata_ %>% filter(sequence_name %in% names(sequence_))
lineB1113_ <- withSeq_ %>% filter(lineage=="B.1.13", !is.na(sample_date))
SeqB.1.13_ <- sequence_[lineB1113_$sequence_name]
#---------------------------------------------------------------- end

#---- compute the dna distance -------------------------------------
#--- //compute the dna distance using SNP cutoff
#--- //extract cluster using cutoff "thrsh_"
#---------------------------------------------------------------- //
distB.1.13_ <- dist.dna(SeqB.1.13_, model="N")
hist(distB.1.13_, xlab="dSNP", main="HIstogram of lineage B.1.13 DNA distance")

clustThrsh <- getClusters(dist=distB.1.13_, thrsh=thrsh_)
clustThrsh1 <- getClusters(dist=distB.1.13_, thrsh=1)
clustThrsh2 <- getClusters(dist=distB.1.13_, thrsh=2)
clustThrsh3 <- getClusters(dist=distB.1.13_, thrsh=3)
#---------------------------------------------------------------- end

#--- //get metadata
metaB.1.13 <- lineB1113_
#--- //set the working directory back to project
setwd(paste(rootProject, "Pairwise Survival Analysis/[Pairwise Survival] COG UK analysis", sep="/"))
#--- //remove unnecessary variables
rm(list=grep("_", ls(), value=T))
rm(rootProject, getClusters)

#---- outputs ---------------------------------------------------
cat("@import [list]'clustThrsh': clusters with cutoff=0 dSNP", "\n")
cat("@import [list]'clustThrsh1': clusters with cutoff=1 dSNP", "\n")
cat("@import [list]'clustThrsh2': clusters with cutoff=2 dSNP", "\n")
cat("@import [list]'clustThrsh3': clusters with cutoff=3 dSNP", "\n")
cat("@import [data.frame]'metaB.1.13': metadata of lineage B.1.13", "\n")
#---------------------------------------------------------------- end
