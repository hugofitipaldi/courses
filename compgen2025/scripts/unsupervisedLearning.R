# Importing the required libraries
library(pheatmap)
library(cluster)
library(ggplot2)
library(viridis)
library(gridExtra)
library(grid)

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

# Calculate the distance matrix between patients using Euclidean distance (default)
d=dist(df)

# Perform hierarchical clustering using complete linkage method
# Complete linkage defines cluster distance as the maximum distance between any two points in different clusters
hc=hclust(d,method="complete")
# Plot the resulting dendrogram showing the hierarchical relationship between patients
plot(hc)

# Load the leukemia gene expression dataset from the compGenomRData package
expFile=system.file("extdata","leukemiaExpressionSubset.rds",package="compGenomRData")
mat=readRDS(expFile) # compGenomRData was manually downloaded from their github page

# Create an annotation data frame for samples
# Extract first 3 characters from column names to identify leukemia types (e.g., "ALL" or "AML")
annotation_col = data.frame(
  LeukemiaType =substr(colnames(mat),1,3))
rownames(annotation_col)=colnames(mat)

# Create a heatmap visualization of the gene expression data
# pheatmap: Pretty Heatmap function with enhanced clustering and annotation capabilities
pheatmap(mat,show_rownames=FALSE,show_colnames=FALSE,annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D2",clustering_distance_cols="euclidean")

# cutting the tree
# Perform hierarchical clustering on the transposed matrix (clustering samples/columns)
hcl=hclust(dist(t(mat)))
# Plot the hierarchical clustering dendrogram
plot(hcl,labels = FALSE, hang= -1)
# Add a red rectangle to dendrogram at height 80 to indicate clusters
rect.hclust(hcl, h = 80, border = "red")
# Cut the dendrogram to create cluster assignments in two different ways:
clu.k5=cutree(hcl,k=5) # Method 1: Cut tree to create exactly 5 clusters
clu.h80=cutree(hcl,h=80) # Method 2: Cut tree at height 80
# Display the number of samples in each cluster when using k=5
table(clu.k5) # number of samples for each cluster

# --- VIDEO 3: k-means clustering and the optimal k ----
# k-means
set.seed(101)

# Perform k-means clustering on transposed matrix (clustering samples/patients, not genes)
# k-means finds k cluster centers and assigns each sample to the nearest center
kclu=kmeans(t(mat),centers=5)  

# Check the number of samples in each cluster
table(kclu$cluster)

# Perform k-medoids clustering (PAM: Partitioning Around Medoids)
# Unlike k-means which uses artificial centroids, k-medoids uses actual data points as cluster centers
kmclu=cluster::pam(t(mat),k=5) #  cluster using k-medoids

# Create data frame combining leukemia type and k-medoids cluster assignment
type2kmclu = data.frame(
  LeukemiaType =substr(colnames(mat),1,3),
  cluster=kmclu$cluster)

table(type2kmclu) # Cross-tabulate leukemia types and cluster assignments for k-medoids

# Create similar data frame for k-means clustering results
type2kclu = data.frame(
  LeukemiaType =substr(colnames(mat),1,3),
  cluster=kclu$cluster)
table(type2kclu) # Cross-tabulate leukemia types and cluster assignments for k-means

# Calculate distance matrix between samples
dists=dist(t(mat))

# Perform multidimensional scaling (MDS) to visualize samples in 2D space
# MDS preserves distances between samples as much as possible
mds=cmdscale(dists)

# plot the patients in the 2D space
plot(mds,pch=19,col=rainbow(5)[kclu$cluster])

# set the legend for cluster colors
legend("bottomright",
       legend=paste("clu",unique(kclu$cluster)),
       fill=rainbow(5)[unique(kclu$cluster)],
       border=NA,box.col=NA)


# Better viz
# Create data frame for plotting
mds_df <- data.frame(
  MDS1 = mds[,1],
  MDS2 = mds[,2],
  Cluster = as.factor(kclu$cluster),
  LeukemiaType = substr(colnames(mat), 1, 3)
)

# Create improved MDS plot
ggplot(mds_df, aes(x = MDS1, y = MDS2, color = Cluster, shape = LeukemiaType)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis_d(option = "D") +
  theme_minimal() +
  labs(title = "MDS Plot of Leukemia Samples",
       subtitle = "Colored by k-means cluster, shaped by leukemia type",
       x = "Dimension 1",
       y = "Dimension 2") +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))


# optimum k
# Determine optimum k using silhouette analysis
set.seed(101)
pamclu=cluster::pam(t(mat),k=5) 
plot(silhouette(pamclu),main=NULL) # Plot silhouette widths for k=5

# Calculate average silhouette width for different values of k (2 through 7)
Ks=sapply(2:7,
          function(i) 
            summary(silhouette(pam(t(mat),k=i)))$avg.width)
# Plot average silhouette width vs. k to identify optimal number of clusters
plot(2:7,Ks,xlab="k",ylab="av. silhouette",type="b",
     pch=19)


# Better Viz
# Function to create silhouette plot for a given k
create_sil_plot <- function(k) {
  pam_result <- pam(t(mat), k = k)
  sil <- silhouette(pam_result)
  sil_df <- data.frame(
    cluster = sil[,1],
    sil_width = sil[,3],
    id = 1:nrow(sil)
  )
  
  # Order by cluster and decreasing silhouette width
  sil_df <- sil_df[order(sil_df$cluster, -sil_df$sil_width),]
  sil_df$id <- 1:nrow(sil_df)
  
  # Calculate average silhouette width
  avg_sil <- mean(sil_df$sil_width)
  
  # Plot
  p <- ggplot(sil_df, aes(x = id, y = sil_width, fill = as.factor(cluster))) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = avg_sil, linetype = "dashed", color = "red") +
    scale_fill_viridis_d(option = "D") +
    labs(title = paste("k =", k),
         subtitle = paste("Average silhouette width:", round(avg_sil, 3)),
         x = "Sample ID (ordered by cluster and silhouette width)",
         y = "Silhouette Width",
         fill = "Cluster") +
    theme_minimal() +
    theme(legend.position = "none") +
    ylim(-0.2, 1)
  
  return(p)
}

# Create plots for k = 2 to 6
sil_plots <- lapply(2:6, create_sil_plot)

# Arrange plots in a grid
grid.arrange(
  grobs = sil_plots,
  ncol = 2,
  top = textGrob("Silhouette Analysis for Different k Values", 
                 gp = gpar(fontface = "bold", fontsize = 14)))

# --- VIDEO 4: Dimensionality reduction ----
#--- PCA on covariance matrix

# Create a subset of the data with only two genes
# Transpose the matrix so samples are in rows and genes in columns
sub.mat=t(mat[rownames(mat) %in% c("ENSG00000100504","ENSG00000105383"),])

# Calculate PCA on this small subset (2 genes across all samples)
# scale() standardizes the data (mean=0, sd=1) before PCA
pr=princomp(scale(sub.mat))

# Visualize the variance explained by each principal component
screeplot(pr)

#----- SVD/PCA (Singular Value Decomposition / Principal Component Analysis)

# Calculate SVD on the scaled expression matrix
d=svd(scale(mat)) 
# Also calculate PCA using prcomp (alternative method)
d1=prcomp(mat,scale. = T)

# Plot the first right singular vector (eigengene/eigenpatient)
plot(d$v[,1],type="b")
# Plot the second right singular vector
plot(d$v[,2],type="b")

# Create heatmap without clustering columns (preserving original sample order)
pheatmap(mat,show_rownames=FALSE,show_colnames=FALSE,cluster_cols = FALSE,
         annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean")






# Project gene expression data onto eigengenes (PCs)
# This transforms data from gene space to PC space
eigengenes=scale(mat) %*% (d$v) # projection on eigenassays
# Visualize first two PCs using smoothScatter (density-based scatter plot)
smoothScatter(eigengenes[,1],eigengenes[,2],pch=19,
              colramp =heat.colors)

# Regular scatter plot of first two PCs
plot(eigengenes[,1],eigengenes[,2],pch=19)

# Project data onto eigenassays (left singular vectors)
assays=t(d$u) %*% scale(mat)
plot(assays[1,],assays[2,],pch=19,
     col=as.factor(annotation_col$LeukemiaType))

