
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ChemmineR")


library("ChemmineR") # Loads the package
library(help="ChemmineR") # Lists all functions and classes 
vignette("ChemmineR") # Opens this PDF manual from R 

install.packages("scatterplot3d")

sdfset <- read.SDFset("PubChem.sdf")

#pdf("Plots.pdf")

hc <- hclust(as.dist(distmat), method="single") 
hc[["labels"]] <- cid(apset) # Assign correct item labels 
plot(as.dendrogram(hc), edgePar=list(col=10, lwd=1), horiz=T)

propma <- atomcountMA(sdfset, addH=FALSE) 
boxplot(propma, col="blue", main="Atom Frequency")


apset <- sdf2ap(sdfset) # Generate atom pair descriptor database for searching
cmp.cluster(db=apset, cutoff = c(0.65, 0.5), quiet=TRUE)[1:4,] # Binning clustering using variable similarity cutoffs. 


cmp.duplicated(apset, type=1)[1:4] # Returns AP duplicates as logical vector 
cmp.duplicated(apset, type=2)[1:4,] # Returns AP duplicates as data frame 

plot(sdfset[c("CMP63","CMP67", "CMP14", "CMP26")], print=FALSE) 

dummy <- cmp.cluster(db=apset, cutoff=0, save.distances="distmat.rda", quiet=TRUE)
load("distmat.rda") 


#clusters <- cmp.cluster(db=apset, cutoff = c(0.7, 0.8, 0.9), quiet = TRUE)
#cluster.visualize(apset, clusters, size.cutoff = 1, quiet = TRUE) # Plots all items.

library(scatterplot3d) 
coord <- cluster.visualize(apset, clusters, size.cutoff=1, dimensions=3, quiet=TRUE) 
scatterplot3d(coord)

#dev.off()****
