# Importing required libraries
library(caret)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RANN)

# --- VIDEO 2: Data pre-processing ----
# Get file paths for data stored in the compGenomRData package
fileLGGexp=system.file("extdata",
                       "LGGrnaseq.rds",
                       package="compGenomRData")
fileLGGann=system.file("extdata",
                       "patient2LGGsubtypes.rds",
                       package="compGenomRData")

# Load gene expression values
# LGG: Lower Grade Glioma RNA-seq data
gexp=readRDS(fileLGGexp)
head(gexp[,1:5])
dim(gexp)

# Load patient annotation data (clinical/subtype information)
patient=readRDS(fileLGGann)
head(patient)

#------- data transformation

par(mfrow=c(1,2))
# Plot histogram of raw gene expression values for the 5th sample
hist(gexp[,5],xlab="gene expression",main="",border="blue4",
     col="cornflowerblue")
# Plot histogram of log-transformed gene expression values
# Adding 1 before log transformation handles zeros (log(0) is undefined)
hist(log10(gexp+1)[,5], xlab="gene expression log scale",main="",
     border="blue4",col="cornflowerblue")

# Apply log transformation to the entire expression matrix
gexp=log10(gexp+1)

# transpose the matrix
tgexp <- t(gexp)

# Better viz
raw_data <- data.frame(
  Expression = gexp[,5],  # Using data before log transformation
  Type = "Raw"
)

log_data <- data.frame(
  Expression = log10(raw_data$Expression + 1),
  Type = "Log10"
)

combined_data <- rbind(raw_data, log_data)

# Create density plots with histograms
ggplot(combined_data, aes(x = Expression, fill = Type)) +
  geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.7, position = "identity") +
  geom_density(alpha = 0.3) +
  facet_wrap(~Type, scales = "free") +
  scale_fill_manual(values = c("Raw" = "#1B9E77", "Log10" = "#D95F02")) +
  theme_minimal() +
  labs(title = "Effect of Log Transformation on Gene Expression Distribution",
       subtitle = "Sample 5",
       x = "Expression Value",
       y = "Density") +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold"))

#----------- data filtering and scaling
# Remove near zero variation variables (genes)
# uniqueCut=15 means if 93.3% (1-1/15) of the values are the same, the variable is removed
# This step helps eliminate genes with very low expression variability
nzv=preProcess(tgexp,method="nzv",uniqueCut = 15)

# Apply the near-zero variance filter using the "predict" function
# Returns the filtered dataset with low-variability genes removed
nzv_tgexp=predict(nzv,tgexp)

# Identify most variable genes based on standard deviation
# Standard deviation is calculated for each column (gene)
SDs=apply(tgexp,2,sd )
# Select the top 1000 genes with highest standard deviation
# order() returns indices sorted by SD values (decreasing=TRUE for highest first)
topPreds=order(SDs,decreasing = TRUE)[1:1000]
# Subset the expression matrix to keep only the top 1000 most variable genes
tgexp=tgexp[,topPreds]

# Center the data (subtract mean from each column)
# This standardizes genes to have a mean of zero
processCenter=preProcess(tgexp, method = c("center"))
tgexp=predict(processCenter,tgexp)

# Create a filter for removing highly correlated variables
# If two genes are correlated above the cutoff (0.9), one is removed
# This reduces redundancy in the dataset
corrFilt=preProcess(tgexp, method = "corr",cutoff = 0.9)
tgexp=predict(corrFilt,tgexp) # Apply the correlation filter

# ---------- dealing with missing values
# Create a copy of the expression data with an artificially introduced NA value
missing_tgexp=tgexp
missing_tgexp[1,1]=NA

anyNA(missing_tgexp) # check if there are NA values

# Method 1: Drop columns (genes) with any missing values
# colSums(is.na()) counts missing values in each column
# Only keep columns where this count is 0 (no missing values)
gexpnoNA=missing_tgexp[ , colSums(is.na(missing_tgexp)) == 0]

# Method 2: Median imputation
# Replace missing values with the median value of that column (gene)
mImpute=preProcess(missing_tgexp,method="medianImpute")
imputedGexp=predict(mImpute,missing_tgexp)

# Method 3: k-nearest neighbors imputation
# Replace missing values based on values from similar samples/genes
knnImpute=preProcess(missing_tgexp,method="knnImpute")
knnimputedGexp=predict(knnImpute,missing_tgexp) 

# --- VIDEO 3: Data split strategies ----
# --------- splitting the data

tgexp=merge(patient,tgexp,by="row.names")

# push sample ids back to the row names
rownames(tgexp)=tgexp[,1]
tgexp=tgexp[,-1]
dim(tgexp)

set.seed(3031) # set the random number seed for reproducibility 

# get indices for 70% of the data set
intrain <- createDataPartition(y = tgexp[,1], p= 0.7)[[1]]

# seperate test and training sets
training <- tgexp[intrain,]
testing <- tgexp[-intrain,]



# -------- predict subtypes with KNN
library(caret)
knnFit=knn3(x=training[,-1], # training set
            y=training[,1], # training set class labels
            k=5)
# predictions on the test set
trainPred=predict(knnFit,training[,-1])
testPred=predict(knnFit,testiing[,-1])


# --------  model performance

# get k-NN prediction on the training data itself, with k=5
knnFit=knn3(x=training[,-1], # training set
            y=training[,1], # training set class labels
            k=5)

# predictions on the test set, return class labels
testPred=predict(knnFit,testing[,-1],type="class")

# compare the predicted labels to real labels
# get different performance metrics
confusionMatrix(data=testing[,1],reference=testPred)



# -------- ROC curves

library(pROC)

# get k-NN class probabilities
# prediction probabilities on the test set
testProbs=predict(knnFit,testing[,-1])

# get the roc curve
rocCurve <- pROC::roc(response = testing[,1],
                      predictor = testProbs[,1],
                      ## This function assumes that the second
                      ## class is the class of interest, so we
                      ## reverse the labels.
                      levels = rev(levels(testing[,1])))
# plot the curve
plot(rocCurve, legacy.axes = TRUE)

# return area under the curve
pROC::auc(rocCurve)



# -----------Random Forests 


set.seed(17)

# we will do no resampling based prediction error
# although it is advised to do so even for random forests
trctrl <- trainControl(method = "none")

# we will now train random forest model
rfFit <- train(subtype~., 
               data = training, 
               method = "ranger",
               trControl=trctrl,
               importance="permutation", # calculate importance
               tuneGrid = data.frame(mtry=100,
                                     min.node.size = 1,
                                     splitrule="gini")
)
# print OOB error
rfFit$finalModel$prediction.error


# ---------- Variable Importance

plot(varImp(rfFit),top=10)


# ---------- Regression with random forests

# file path for CpG methylation and age
fileMethAge=system.file("extdata",
                        "CpGmeth2Age.rds",
                        package="compGenomRData")

# read methylation-age table
ameth=readRDS(fileMethAge)
dim(ameth)

# filter based on variance
ameth=ameth[,c(TRUE,matrixStats::colSds(as.matrix(ameth[,-1]))>0.1)]
dim(ameth)

# we are not going to do any cross-validatin
# and rely on OOB error
trctrl <- trainControl(method = "none")

# we will now train random forest model
rfregFit <- train(Age~., 
                  data = ameth, 
                  method = "ranger",
                  trControl=trctrl,
                  # calculate importance
                  importance="permutation", 
                  tuneGrid = data.frame(mtry=50,
                                        min.node.size = 5,
                                        splitrule="variance")
)
# plot Observed vs OOB predicted values from the model
plot(ameth$Age,rfregFit$finalModel$predictions,
     pch=19,xlab="observed Age",
     ylab="OOB predicted Age")
mtext(paste("R-squared",
            format(rfregFit$finalModel$r.squared,digits=2)))

