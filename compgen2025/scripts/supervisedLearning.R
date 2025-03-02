# Importing required libraries
library(caret)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RANN)
library(pROC)
library(ggplot2)
library(viridis)
library(dplyr)

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
# Merge patient annotation data with gene expression data
# by="row.names" joins the datasets on their row names (patient IDs)
tgexp=merge(patient,tgexp,by="row.names")

# Push sample IDs back to the row names
# After merge, the row names become Row.names column, so we restore them
rownames(tgexp)=tgexp[,1]
tgexp=tgexp[,-1]
dim(tgexp)

set.seed(3031) # set the random number seed for reproducibility 

# Create indices for 70% of the data for training
# createDataPartition ensures balanced sampling across classes in column 1
intrain <- createDataPartition(y = tgexp[,1], p= 0.7)[[1]]

# Split data into training (70%) and testing (30%) sets
training <- tgexp[intrain,]
testing <- tgexp[-intrain,]

# --- VIDEO 4: Training and model assessment ----
# -------- predict subtypes with KNN

# Train KNN model with k=5
knnFit = knn3(x = training[,-1],     # Training features (gene expression data)
              y = training[,1],      # Training class labels (subtypes)
              k = 5)                 # Number of neighbors to consider

# Make predictions on both training and test sets
trainPred=predict(knnFit,training[,-1])
testPred=predict(knnFit,testing[,-1])

# --------  model performance

# Train a k-NN classifier on the training data
knnFit=knn3(x=training[,-1], # Training features (gene expression data)
            y=training[,1], # Training class labels (subtypes)
            k=5) # Number of neighbors to consider

# Make predictions on the test set
# type="class" returns predicted class labels directly instead of probabilities
testPred=predict(knnFit,testing[,-1],type="class")

# Evaluate model performance by comparing predicted vs. actual labels
# confusionMatrix calculates various performance metrics
confusionMatrix(data=testing[,1],reference=testPred)

# -------- ROC curves

# Get k-NN prediction probabilities for the test set
testProbs=predict(knnFit,testing[,-1])

# Calculate the ROC curve
rocCurve <- pROC::roc(response = testing[,1],        # Actual class labels
                      predictor = testProbs[,1],     # Predicted probabilities for the first class
                      # This function assumes that the second class is the class of interest,
                      # so we reverse the levels to make the first class the positive class
                      levels = rev(levels(testing[,1])))

# plot the curve
par(mfrow = c(1, 1))
plot(rocCurve, legacy.axes = TRUE)

# return area under the curve
pROC::auc(rocCurve)

# --- VIDEO 6: Random Forests ----
# -----------Random Forests 

set.seed(17)

# Define training control parameters
# method="none" means no cross-validation or resampling will be used
# Although resampling is generally recommended even for random forests
trctrl <- trainControl(method = "none")

# Train a Random Forest model using the ranger implementation
rfFit <- train(subtype~.,           # Formula: predict subtype using all other variables
               data = training,     # Training dataset
               method = "ranger",   # Use ranger implementation of Random Forests
               trControl = trctrl,  # Training control parameters defined above
               importance = "permutation",  # Calculate variable importance using permutation
               tuneGrid = data.frame(mtry = 100,           # Number of variables to consider at each split
                                     min.node.size = 1,     # Minimum node size
                                     splitrule = "gini")    # Split rule criterion
)

# Print Out-Of-Bag (OOB) error rate
# This is an unbiased estimate of the test error
rfFit$finalModel$prediction.error

# --- VIDEO 7: Variable importance ----
# ---------- Variable Importance

plot(varImp(rfFit),top=10)

# Better Viz
# Extract variable importance
var_imp <- as.data.frame(ranger::importance(rfFit$finalModel))
var_imp$Variable <- rownames(var_imp)
names(var_imp)[1] <- "Importance"

# Sort by importance and get top 20
top_vars <- var_imp %>%
  arrange(desc(Importance)) %>%
  head(20)

# Create importance plot
ggplot(top_vars, aes(x = reorder(Variable, Importance), y = Importance, fill = Importance)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 Variables by Importance in Random Forest",
       subtitle = "Based on permutation importance",
       x = "",
       y = "Variable Importance") +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))

# --- VIDEO 8: Regression with ML ----
# ---------- Regression with random forests

# Get file path for CpG methylation and age data
fileMethAge=system.file("extdata",
                        "CpGmeth2Age.rds",
                        package="compGenomRData")

# Read methylation-age table
ameth=readRDS(fileMethAge)
dim(ameth)

# Filter CpG sites based on variance
# Keep Age column (first column) and CpGs with standard deviation > 0.1
ameth=ameth[,c(TRUE,matrixStats::colSds(as.matrix(ameth[,-1]))>0.1)]
dim(ameth) # Display dimensions after filtering

# Set up training control
# We're not doing cross-validation and will rely on OOB error instead
trctrl <- trainControl(method = "none")

# Train random forest regression model
rfregFit <- train(Age ~ .,             # Predict Age using all other variables
                  data = ameth,         # Methylation dataset
                  method = "ranger",    # Use ranger implementation of Random Forests
                  trControl = trctrl,   # Training control parameters
                  # Calculate variable importance using permutation method
                  importance = "permutation",
                  # Set hyperparameters
                  tuneGrid = data.frame(mtry = 50,              # Number of variables to consider at each split
                                        min.node.size = 5,       # Minimum node size
                                        splitrule = "variance"))  # Split rule criterion for regression


# Plot Observed vs. OOB predicted values from the model
plot(ameth$Age,rfregFit$finalModel$predictions,
     pch=19,xlab="observed Age",
     ylab="OOB predicted Age")
mtext(paste("R-squared",
            format(rfregFit$finalModel$r.squared,digits=2))) # Add R-squared value to the plot

# Better Viz
# Create data frame for plotting
pred_df <- data.frame(
  Observed = ameth$Age,
  Predicted = rfregFit$finalModel$predictions,
  Error = rfregFit$finalModel$predictions - ameth$Age
)

# Calculate statistics
rmse <- sqrt(mean((pred_df$Observed - pred_df$Predicted)^2))
r_squared <- rfregFit$finalModel$r.squared
mae <- mean(abs(pred_df$Error))

# Create enhanced scatter plot
ggplot(pred_df, aes(x = Observed, y = Predicted)) +
  geom_point(aes(color = abs(Error)), alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "solid") +
  scale_color_viridis(option = "D") +
  theme_minimal() +
  labs(title = "Random Forest Regression: Predicted vs. Observed Age",
       subtitle = paste("R² =", round(r_squared, 3), "| RMSE =", round(rmse, 2), "| MAE =", round(mae, 2)),
       x = "Observed Age (years)",
       y = "Predicted Age (years)",
       color = "Absolute Error") +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

# Add additional analysis to prediction data frame
pred_df$Residual <- pred_df$Error
pred_df$StdResidual <- scale(pred_df$Residual)

# Create residual plots
p1 <- ggplot(pred_df, aes(x = Predicted, y = Residual)) +
  geom_point(aes(color = abs(Residual)), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Residual Plot",
       x = "Predicted Age",
       y = "Residual (Predicted - Observed)",
       color = "Absolute\nResidual") +
  theme(legend.position = "none")

p2 <- ggplot(pred_df, aes(x = Residual)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Residual Distribution",
       x = "Residual (Predicted - Observed)",
       y = "Density")

# Arrange plots
gridExtra::grid.arrange(p1, p2, ncol = 2)
