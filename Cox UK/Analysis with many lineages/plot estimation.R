#---- //Plots ----
#     * The estimation of alpha, beta, and mixture mean
#     * Plot the distribution of mixture mean from bootstrapping
#----------------------------------------------------------------------//
require(ggplot2)
alpha_ <- ggplot(record, aes(x=cutoff)) + 
  geom_boxplot(aes(y=alpha, fill=ncomp), alpha=.3) +
  theme(legend.position="bottom") + xlab("") +
  ggtitle("Lineage B")
beta_ <- ggplot(record, aes(x=cutoff)) + 
  geom_boxplot(aes(y=beta, fill=ncomp), alpha=.3) +
  theme(legend.position="bottom") + xlab("") +
  ggtitle("Lineage B")
mixtureMean_ <- ggplot(record, aes(x=cutoff)) + 
  geom_boxplot(aes(y=mixtureMean, fill=ncomp), alpha=.3) +
  theme(legend.position="bottom") + xlab("") +
  ggtitle("Lineage B")

#--- //plot the dist. of mean of serial interval with ncomp=4
hist_ <- ggplot(record %>% filter(ncomp=="ncomp=4")) +
  geom_histogram(aes(x=mixtureMean, y=..density..), fill="red", alpha=.3, col="black", bins=30) +
  ggtitle("Dist. of mean of serial interval from bootstrapping", subtitle="ncomp=4") +
  facet_wrap(~cutoff) + xlab("")
