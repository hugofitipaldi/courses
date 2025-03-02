# Loading packages
library(ggplot2)
library(dplyr)
library(mosaic)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("qvalue")
library(qvalue)
require(matrixStats)
library(gridExtra)
library(tidyr)  

# --- VIDEO 1: Summarizing data ----
# --- Central Tendency and spread ---

x=rnorm(20,mean=6,sd=0.7)

mean(x)
median(x)

var(x)
sd(x)
IQR(x)

# Creating the plot for the next exercise
# Set parameters for the normal distribution
mu <- 0
sigma <- 2

# Create a data frame for the curve
x <- seq(-6, 6, length.out = 1000)
y <- dnorm(x, mean = mu, sd = sigma)
df <- data.frame(x = x, y = y)

# Calculate the specific x value (-2) probability area
x_cut <- -2
prob_area <- pnorm(x_cut, mean = mu, sd = sigma)

# Plot
ggplot(df, aes(x = x, y = y)) +
  geom_line() +
  geom_area(data = df, aes(x = x, y = y), fill = "#a3d977", alpha = 0.7) +
  geom_area(data = df %>% filter(x <= x_cut), 
            aes(x = x, y = y), fill = "#432371", alpha = 0.8) +
  geom_vline(xintercept = -2, linetype = "solid", color = "black") +
  annotate("text", x = -1, y = 0.21, label = "z = -1", vjust = 1.5) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "#e5e5e5"),
    panel.grid.minor = element_line(color = "#e5e5e5"),
    plot.background = element_rect(fill = "#f5f5f5", color = NA)) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-6, 7)) +
  scale_y_continuous(limits = c(0, 0.22)) +
  labs(x = "x", y = "density")

#get the probability of P(X =< -2) | mean=0 and sd=2 
pnorm(-2, mean=0, sd=2)

#get the probability of P(X > -2) | mean=0 and sd=2 
pnorm(-2, mean=0, sd=2,lower.tail = FALSE)

#get 5 random numbers from norm. dist. | mean=0 and sd=2
rnorm(5, mean=0 , sd=2)

# --- VIDEO 2: Precision of estimates: Confidence Intervals ----
#--- precision of estimates ---

set.seed(21)
sample1= rnorm(50,mean=20,sd=5) # simulate a sample

# do bootstrap resampling, sampling with replacement
boot.means=do(1000) * mean(resample(sample1))

# get percentiles from the bootstrap means
q=quantile(boot.means[,1],p=c(0.025,0.975))

# plot the histogram
hist(boot.means[,1],col="cornflowerblue",border="white",
     xlab="sample means")
abline(v=c(q[1], q[2] ),col="red")
text(x=q[1],y=200,round(q[1],3),adj=c(1,0))
text(x=q[2],y=200,round(q[2],3),adj=c(0,0))

# a better version of this plot:
boot_df <- data.frame(means = boot.means[,1])

ggplot(boot_df, aes(x = means)) +
  geom_histogram(
    fill = "cornflowerblue", 
    color = "white",
    alpha = 0.7,
    bins = 30) +
  geom_vline(
    xintercept = q[1], 
    color = "#FF5555", 
    size = 1.2,
    linetype = "dashed") +
  geom_vline(
    xintercept = q[2], 
    color = "#FF5555", 
    size = 1.2,
    linetype = "dashed") +
  annotate(
    "text",
    x = q[1] - 0.1,
    y = 100,
    label = round(q[1], 3),
    color = "#FF5555",
    fontface = "bold",
    hjust = 1) +
  annotate(
    "text",
    x = q[2] + 0.1,
    y = 100,
    label = round(q[2], 3),
    color = "#FF5555",
    fontface = "bold",
    hjust = 0) +
  labs(
    title = "Bootstrap Distribution of Sample Means",
    subtitle = paste0("95% CI: [", round(q[1], 3), ", ", round(q[2], 3), "]"),
    x = "Sample Means",
    y = "Frequency") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "darkgrey"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"))

alpha=0.05 # 95% confidence interval
sd=5
n=50
mean(sample1)+qnorm(c(alpha/2,1-alpha/2))*sd/sqrt(n) # Similar to the one calculated in the bootstrap

# --- VIDEO 3: Testing differences between samples: Hypothesis testing ----
#--- hypothesis testing ---

set.seed(100)
gene1=rnorm(30,mean=4,sd=2)
gene2=rnorm(30,mean=2,sd=2)
org.diff=mean(gene1)-mean(gene2)
gene.df=data.frame(exp=c(gene1,gene2),
                   group=c( rep("test",30),rep("control",30) ) )


exp.null <- do(1000) * diff(mosaic::mean(exp ~ shuffle(group), data=gene.df))

hist(exp.null[,1],xlab="null distribution | no difference in samples",
     main=expression(paste(H[0]," :no difference in means") ),
     xlim=c(-2,2),col="cornflowerblue",border="white")
abline(v=quantile(exp.null[,1],0.95),col="red" )
abline(v=org.diff,col="blue" )
text(x=quantile(exp.null[,1],0.95),y=200,"0.05",adj=c(1,0),col="red")
text(x=org.diff,y=200,"org. diff.",adj=c(1,0),col="blue")

p.val=sum(exp.null[,1]>org.diff)/length(exp.null[,1])
p.val

# Better plot:
null_dist_df <- data.frame(diff_means = exp.null[,1])
crit_val <- quantile(exp.null[,1], 0.95)

# Create the plot
ggplot(null_dist_df, aes(x = diff_means)) +
  geom_histogram(
    fill = "cornflowerblue", 
    color = "white",
    alpha = 0.8,
    bins = 25) +
  geom_vline(
    xintercept = crit_val, 
    color = "#FF5555", 
    size = 1,
    linetype = "dashed") +
  geom_vline(
    xintercept = org.diff, 
    color = "#4477AA", 
    size = 1,
    linetype = "solid") +
  annotate(
    "text",
    x = crit_val - 0.05,
    y = 100,
    label = "0.05",
    color = "#FF5555",
    fontface = "bold",
    hjust = 1) +
  annotate(
    "text",
    x = org.diff - 0.05,
    y = 125,
    label = "org. diff.",
    color = "#4477AA",
    fontface = "bold",
    hjust = 1) +
  annotate(
    "text",
    x = max(c(org.diff, crit_val)) + 0.3,
    y = 75,
    label = paste("p-value =", round(p.val, 4)),
    color = "black",
    fontface = "italic",
    hjust = 0) +
  geom_area(
    data = subset(null_dist_df, diff_means > crit_val),
    aes(x = diff_means, y = after_stat(count)),
    stat = "bin",
    bins = 25,
    fill = "#FF5555",
    alpha = 0.3) +
  labs(
    title = expression(paste(H[0], ": No difference in means")),
    subtitle = "Permutation-based null distribution",
    x = "Difference in means (null distribution)",
    y = "Frequency",
    caption = "Red area: rejection region (α = 0.05)") +
  xlim(-2, 2) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "darkgrey"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"))

# Welch's t-test (default in R)
# This version does NOT assume equal variances between groups
# It adjusts the degrees of freedom using the Welch-Satterthwaite equation
# This is generally recommended as it's more robust when variances differ
stats::t.test(gene1, gene2) # default in R for t-test

# Student's t-test with equal variance assumption
# This is the classical t-test that assumes both groups have equal variances
# Only use this when you have strong reason to believe variances are equal
# or when sample sizes are very small and approximately equal
stats::t.test(gene1, gene2, var.equal=TRUE)

# --- VIDEO 4: multiple testing correction ----
#--- mult. test correction ---

data(hedenfalk)

qvalues = qvalue(hedenfalk$p)$q
bonf.pval=p.adjust(hedenfalk$p,method ="bonferroni")
fdr.adj.pval=p.adjust(hedenfalk$p,method ="fdr")

plot(hedenfalk$p,qvalues,pch=19,ylim=c(0,1),
     xlab="raw P-values",ylab="adjusted P-values")
points(hedenfalk$p,bonf.pval,pch=19,col="red")
points(hedenfalk$p,fdr.adj.pval,pch=19,col="blue")
legend("bottomright",legend=c("q-value","FDR (BH)","Bonferroni"),
       fill=c("black","blue","red"))

# Better plot:
plot_data <- data.frame(
  raw_p = hedenfalk$p,
  qvalue = qvalues,
  bonferroni = bonf.pval,
  fdr = fdr.adj.pval
)

plot_data_long <- tidyr::pivot_longer(
  plot_data, 
  cols = c(qvalue, bonferroni, fdr),
  names_to = "adjustment_method",
  values_to = "adjusted_p"
)

plot_data_long$adjustment_method <- factor(
  plot_data_long$adjustment_method,
  levels = c("qvalue", "fdr", "bonferroni"),
  labels = c("q-value", "FDR (BH)", "Bonferroni")
)

ggplot(plot_data_long, aes(x = raw_p, y = adjusted_p, color = adjustment_method)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_smooth(se = FALSE, linewidth = 0.7, alpha = 0.7, method = "loess") +
  scale_color_manual(
    values = c("q-value" = "#1b9e77", "FDR (BH)" = "#7570b3", "Bonferroni" = "#d95f02"),
    name = "Adjustment Method"
  ) +
  labs(
    title = "Comparison of P-value Adjustment Methods",
    subtitle = "Hedenfalk dataset",
    x = "Raw P-values",
    y = "Adjusted P-values",
    caption = "Dashed line: y=x reference"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = "gray90"),
    legend.margin = margin(5, 5, 5, 5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray90", fill = NA),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
    axis.title = element_text(face = "bold")) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))

# --- VIDEO 5: moderated tests ----
#--- moderated t-test ---
# The moderated t-test is a variation of the standard t-test that "shrinks" individual gene variances 
# toward a common value to improve variance estimates when sample sizes are small (common in genomics)

set.seed(100)

# Generate sample data matrix from normal distribution
# Creating 3000 values with mean=200 and sd=70 to simulate gene expression levels
gset = rnorm(3000, mean=200, sd=70)
data = matrix(gset, ncol=6)  # Reshape into matrix with 6 columns (samples) and 500 rows (genes)

# Define sample groups
group1 = 1:3  # First 3 columns represent group 1 (e.g., treatment)
group2 = 4:6  # Last 3 columns represent group 2 (e.g., control)
n1 = 3  # Sample size of group 1
n2 = 3  # Sample size of group 2

# Calculate mean difference between groups for each gene
dx = rowMeans(data[,group1]) - rowMeans(data[,group2])

# Calculate the standard error using pooled variance estimate
# This is the denominator of the t-statistic using classical approach
stderr <- sqrt((rowVars(data[,group1])*(n1-1) + rowVars(data[,group2])*(n2-1)) / 
                 (n1+n2-2) * (1/n1 + 1/n2))

# Moderate (shrink) the standard error by averaging with the median standard error
# This is the key step in the moderated t-test approach that stabilizes variance estimates
# Particularly useful when sample sizes are small, as in many genomics experiments
mod.stderr <- (stderr + median(stderr)) / 2  # moderation in variation

# Calculate t-statistic using moderated standard error
t.mod = dx / mod.stderr

# Calculate two-sided p-values for moderated t-test
# pt() gives the cumulative distribution function for t-distribution
# We multiply by 2 for a two-sided test
p.mod = 2*pt(-abs(t.mod), n1+n2-2)

# Calculate regular t-statistic using unmoderated standard error
t = dx / stderr

# Calculate two-sided p-values for regular t-test
p = 2*pt(-abs(t), n1+n2-2)

# Create side-by-side histograms to compare results
par(mfrow=c(1,2))  # Set up plotting area for 2 plots side by side

# Plot histogram of regular t-test p-values
hist(p, col="cornflowerblue", border="white", main="", xlab="P-values t-test")
mtext(paste("significant tests:", sum(p<0.05)))  # Add text for number of significant results

# Plot histogram of moderated t-test p-values
hist(p.mod, col="cornflowerblue", border="white", main="", xlab="P-values mod. t-test")
mtext(paste("significant tests:", sum(p.mod<0.05)))

# Better Viz
# Reproducing the calculation from above
set.seed(100)
gset = rnorm(3000, mean=200, sd=70)
data = matrix(gset, ncol=6)
group1 = 1:3
group2 = 4:6
n1 = 3
n2 = 3

# Calculating mean difference between groups
dx = rowMeans(data[,group1]) - rowMeans(data[,group2])

# Calculating standard errors
stderr <- sqrt((rowVars(data[,group1])*(n1-1) + rowVars(data[,group2])*(n2-1)) / 
                 (n1+n2-2) * (1/n1 + 1/n2))
mod.stderr <- (stderr + median(stderr)) / 2  # moderation in variation

# Calculating t-statistics and p-values
t.mod = dx / mod.stderr
p.mod = 2*pt(-abs(t.mod), n1+n2-2)
t = dx / stderr
p = 2*pt(-abs(t), n1+n2-2)

# df for plotting
results_df <- data.frame(
  GeneID = 1:length(p),
  Standard_P = p,
  Moderated_P = p.mod,
  Standard_t = t,
  Moderated_t = t.mod,
  Standard_Significant = p < 0.05,
  Moderated_Significant = p.mod < 0.05
)

p_values_long <- pivot_longer(
  results_df, 
  cols = c(Standard_P, Moderated_P),
  names_to = "Test_Type", 
  values_to = "P_Value"
)

p_values_long$Test_Type <- factor(
  p_values_long$Test_Type,
  levels = c("Standard_P", "Moderated_P"),
  labels = c("Standard t-test", "Moderated t-test")
)

# Counting significant tests for each method
sig_counts <- p_values_long %>%
  group_by(Test_Type) %>%
  summarize(
    Significant = sum(P_Value < 0.05),
    Total = n(),
    Percentage = round(Significant/Total*100, 1)
  )

# Create plots
# 1. Histogram comparison
p1 <- ggplot(p_values_long, aes(x = P_Value, fill = Test_Type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Standard t-test" = "#4477AA", "Moderated t-test" = "#CC6677")) +
  labs(
    title = "Comparison of P-value Distributions",
    subtitle = paste(
      "Standard t-test:", sig_counts$Significant[1], "significant genes (", sig_counts$Percentage[1], "%)\n",
      "Moderated t-test:", sig_counts$Significant[2], "significant genes (", sig_counts$Percentage[2], "%)"
    ),
    x = "P-value",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# 2. Scatter plot comparing p-values from both methods
p2 <- ggplot(results_df, aes(x = Standard_P, y = Moderated_P)) +
  geom_point(alpha = 0.3, aes(color = Standard_P != Moderated_P)) +
  geom_abline(linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "red") +
  scale_color_manual(values = c("TRUE" = "#CC6677", "FALSE" = "gray70"), guide = "none") +
  labs(
    title = "Standard vs. Moderated t-test P-values",
    x = "Standard t-test P-value",
    y = "Moderated t-test P-value"
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

grid.arrange(p1, p2, ncol = 2)

# Create volcano plots for both methods
# 1. Standard t-test volcano plot
p3 <- ggplot(results_df, aes(x = dx, y = -log10(Standard_P))) +
  geom_point(aes(color = Standard_Significant), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "#CC6677", "FALSE" = "gray70")) +
  labs(
    title = "Standard t-test Volcano Plot",
    x = "Mean Difference",
    y = "-log10(P-value)",
    color = "Significant"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# 2. Moderated t-test volcano plot
p4 <- ggplot(results_df, aes(x = dx, y = -log10(Moderated_P))) +
  geom_point(aes(color = Moderated_Significant), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "#CC6677", "FALSE" = "gray70")) +
  labs(
    title = "Moderated t-test Volcano Plot",
    x = "Mean Difference",
    y = "-log10(P-value)",
    color = "Significant"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

grid.arrange(p1, p2, p3, p4, ncol = 2)

#---- regression ----



# set random number seed, so that the random numbers from the text
# is the same when you run the code.
set.seed(32)

# get 50 X values between 1 and 100
x = runif(50,1,100)

# set b0,b1 and variance (sigma)
b0 = 100
b1 = 2
sigma = 20
# simulate error terms from normal distribution
eps = rnorm(50,0,sigma)
# get y values from the linear equation and addition of error terms
y = b0 + b1*x+ eps

mod1=lm(y~x)

# plot the data points
plot(x,y,pch=20,
     ylab="Gene Expression",xlab="Histone modification score")
# plot the linear fit
abline(mod1,col="blue")

mod1


summary(mod1)

# get confidence intervals 
confint(mod1)

# pull out coefficients from the model
coef(mod1)



