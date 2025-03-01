# Load necessary libraries
library(data.table)
library(stats)

# Read the data
data <- fread("../data/leukemiaExp.txt")

# Prepare the data
gene_names <- data[[1]]  # First column contains gene names
control_data <- as.matrix(data[, 2:31])  # Columns 1-30 are control samples
leukemia_data <- as.matrix(data[, 32:61])  # Columns 31-60 are leukemia samples

# Perform t-test for each gene
t_test_results <- apply(cbind(control_data, leukemia_data), 1, function(x) {
  t.test(x[1:30], x[31:60])
})

# Extract p-values
p_values <- sapply(t_test_results, function(x) x$p.value)

# Adjust p-values for multiple testing
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Identify significantly differentially expressed genes
significant_genes <- gene_names[adjusted_p_values < 0.05]

# Create a data frame with results
results <- data.frame(
  Gene = gene_names,
  P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values,
  Significant = adjusted_p_values < 0.05
)

# Sort results by adjusted p-value
results <- results[order(results$Adjusted_P_Value), ]

# Print the number of significant genes
cat("Number of significantly differentially expressed genes:", sum(results$Significant), "\n")

# Print the top 10 most significant genes
print(head(results, 10))

# Write results to a file
write.csv(results, "../data/differential_expression_results.csv", row.names = FALSE)
