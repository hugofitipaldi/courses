# Importing the required libraries
library(pheatmap)

# --- VIDEO 1: Distance metrics ----
# Create a data frame with gene expression values for 3 genes across 4 patients
df = data.frame(
  IRX4 = c(11, 13, 2, 1),    # Expression values for IRX4 gene
  OCT4 = c(10, 13, 4, 3),    # Expression values for OCT4 gene
  PAX6 = c(1, 3, 10, 9),     # Expression values for PAX6 gene
  row.names = c("patient1", "patient2", "patient3", "patient4")  # Patient identifiers as row names
)
df
# Calculate distances between patients using different metrics:

# Manhattan distance (L1 norm): sum of absolute differences
# This measures the sum of absolute differences in expression values
# Also known as city block or taxicab distance
dist(df, method = "manhattan")

# Euclidean distance (L2 norm): square root of sum of squared differences
# This is the "straight line" distance in multidimensional space
dist(df, method = "euclidean")

# Correlation distance: 1 minus the correlation coefficient
# Measures similarity in expression patterns rather than absolute values
# Transposing with t() because we want correlations between patients (rows), not genes
as.dist(1 - cor(t(df)))

# scaling data 
# The scale() function standardizes each column (variable) in the data frame
# This transforms values to have mean = 0 and standard deviation = 1
scale(df)

# --- VIDEO 2: Hierarchical clustering ----
##  Hierarchical clustering

d=dist(df)
hc=hclust(d,method="complete")
plot(hc)


expFile=system.file("extdata","leukemiaExpressionSubset.rds",package="compGenomRData")
mat=readRDS(expFile)

# set the leukemia type annotation for each sample
annotation_col = data.frame(
  LeukemiaType =substr(colnames(mat),1,3))
rownames(annotation_col)=colnames(mat)


pheatmap(mat,show_rownames=FALSE,show_colnames=FALSE,annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D2",clustering_distance_cols="euclidean")


# cutting the tree

hcl=hclust(dist(t(mat)))
plot(hcl,labels = FALSE, hang= -1)
rect.hclust(hcl, h = 80, border = "red")

clu.k5=cutree(hcl,k=5) # cut tree so that there are 5 clusters
clu.h80=cutree(hcl,h=80) # cut tree/dendrogram from height 80

table(clu.k5) # number of samples for each cluster





# k-means
set.seed(101)

# we have to transpose the matrix t()
# so that we calculate distances between patients
kclu=kmeans(t(mat),centers=5)  

# number of data points in each cluster
table(kclu$cluster)

kmclu=cluster::pam(t(mat),k=5) #  cluster using k-medoids

# make a data frame with Leukemia type and cluster id
type2kmclu = data.frame(
  LeukemiaType =substr(colnames(mat),1,3),
  cluster=kmclu$cluster)

table(type2kmclu)


type2kclu = data.frame(
  LeukemiaType =substr(colnames(mat),1,3),
  cluster=kclu$cluster)
table(type2kclu)

# Calculate distances
dists=dist(t(mat))

# calculate MDS
mds=cmdscale(dists)

# plot the patients in the 2D space
plot(mds,pch=19,col=rainbow(5)[kclu$cluster])

# set the legend for cluster colors
legend("bottomright",
       legend=paste("clu",unique(kclu$cluster)),
       fill=rainbow(5)[unique(kclu$cluster)],
       border=NA,box.col=NA)



# optimum k

library(cluster)
set.seed(101)
pamclu=cluster::pam(t(mat),k=5)
plot(silhouette(pamclu),main=NULL)

#plot silhouette for different k
Ks=sapply(2:7,
          function(i) 
            summary(silhouette(pam(t(mat),k=i)))$avg.width)
plot(2:7,Ks,xlab="k",ylab="av. silhouette",type="b",
     pch=19)


#--- PCA on covariance matrix

# create the subset of the data with two genes only
# notice that we transpose the matrix so samples are 
# on the columns
sub.mat=t(mat[rownames(mat) %in% c("ENSG00000100504","ENSG00000105383"),])


# calculate the PCA only for our genes and all the samples
pr=princomp(scale(sub.mat))

# show explained variance
screeplot(pr)


#----- SVD/PCA

d=svd(scale(mat)) 
d1=prcomp(mat,scale. = T)

# plot first eigengene/vector
plot(d$v[,1],type="b")
plot(d$v[,2],type="b")

# plot the heatmap without clustering columns
pheatmap(mat,show_rownames=FALSE,show_colnames=FALSE,cluster_cols = FALSE,
         annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean")






# projections on eigengenes or assays
eigengenes=scale(mat) %*% (d$v) # projection on eigenassays
smoothScatter(eigengenes[,1],eigengenes[,2],pch=19,
              colramp =heat.colors)

plot(eigengenes[,1],eigengenes[,2],pch=19)

assays=t(d$u) %*% scale(mat) # projection on eigenassays
plot(assays[1,],assays[2,],pch=19,
     col=as.factor(annotation_col$LeukemiaType))

