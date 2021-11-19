# PCA Code from Abrahms et al. 2017 Movement Ecology
# Based on tutorial located at 
# https://tgmstat.wordpress.com/2013/11/28/computing-and-visualizing-pca-in-r/#ref1
source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)


# Data should be in dataframe or matrix format
Data<-PCA_Data[,c("VI", "T2R", "MNSD", "RT", "TAC")] 
Data<-log(Data) # log transform all data
spp<-PCA_Data[,"SpeciesName"]
Classification<-PCA_Data[,"Classification"]
Marker<-PCA_Data[,"Marker"]

# Use function prcomp
# center = TRUE (default setting) centers the data so that each column has a mean of 0
# scale = TRUE (default is false) scales the data so each column has a stdev of 1
# tol=n gives the PCs whose standard deviations are greater than n% of the standard deviation of the first PC
# deviation of the first principal component
pca.full=prcomp(Data, scale=T)
barplot(pca.full$sdev/pca.full$sdev[1]) #visualize stdev cutoff for tol parameter
pca=prcomp(Data, scale=T, tol=0.6)
pca
summary(pca) # use summary to see cumulative contribution of each PC to the variance

# Plot eigenvalues and percentages of variation of an ordination object
# from: http://rfunctions.blogspot.com/2015/01/pca-principal-component-analysis.html
# Kaiser rule and broken stick model
# Usage:
# evplot(ev)
# where ev is a vector of eigenvalues
# Author: Francois Gillet, 25 August 2012
ev <- pca.full$sdev^2
evplot <- function(ev)
{
  # Broken stick model (MacArthur 1957)
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}
evplot(ev)


# Plot data on PCA axes
g <- ggbiplot(pca, choices=1:2, scale=1, obs.scale = 1, var.scale = 1, groups = spp, ellipse = F, ellipse.prob = 0.5)
g <- g + scale_shape_manual(values=c(17,4,16)) + geom_point(aes(shape = Marker, colour=spp), size=2.5) + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'vertical', legend.position = 'right')
g <- g + ylim(-3,3)
g <- g + theme_bw() + theme(panel.grid=element_blank()) + theme(legend.title=element_blank())
plot(g)

p <- ggbiplot(pca, choices=1:2, scale=1, obs.scale = 1, var.scale = 1, groups = Classification, ellipse = T, ellipse.prob = 0.5)
p <- p + scale_shape_manual(values=c(17,4,16)) + geom_point(aes(shape = Marker, colour=Classification), size=2.5) + scale_color_discrete(name = '')
p <- p + theme(legend.direction = 'vertical', legend.position = 'right')
p <- p + scale_colour_brewer(palette = "Spectral")
p <- p + ylim(-3,3)
p <- p + theme_bw() + theme(panel.grid=element_blank()) + theme(legend.title=element_blank())
plot(p)


# CLUSTER ANALYSIS
# Determine number of clusters - look for bend in graph
# from B.S. Everitt & T. Hothorn, A Handbook of Statistical Analyses Using R (2nd ed.) pg 251
wss <- (nrow(Data)-1)*sum(apply(Data,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(Data, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

# Ward Hierarchical Cluster Analysis
D <- pca.full$x[,1:3] # Use PCA values for cluster analysis
row.names(D)<-PCA_Data$ID
d <- dist(D, method = "euclidean") # distance matrix
dend <- hclust(d, method="ward.D")
A2Rplot(dend, k = 4, boxes = FALSE, col.up = "gray50",lty.up=1, col.down = c("#FDAE61","#D7191C","#ABDDA4", "#2B83BA"), main=NULL, cex.axis=0.1)

# Permutation test (bootstr) to generate cluster p-values
library(pvclust)
fit <- pvclust(D, method.hclust="ward.D", method.dist="euclidean")
plot(fit)
pvrect(fit, alpha=.95)
print(fit, digits=3)
msplot(fit, edges=c(147,148,149,150))




