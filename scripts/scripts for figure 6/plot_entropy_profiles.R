# Required packages:
# ggplot2

library(ggplot2)
library(readr)

#############################################################
# Shannon entropy plots for miR-548 family sequence analysis
#############################################################

# Load data
data_3p <- read.csv("3p-sequence_analysis_results.csv")

# Plot for miR-548-3p
ggplot(data_3p, aes(x = Position, y = Entropy)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = c(5,11), linetype="dashed") +
  annotate("rect", xmin=5, xmax=11, ymin=-Inf, ymax=Inf, alpha=0.2) +
  theme_minimal() +
  labs(title = "Shannon Entropy per Position (miR-548-3p sequences)",
       x = "Position",
       y = "Entropy")


# Load data
data_5p <- read.csv("5p-sequence_analysis_results.csv")

# Plot for miR-548-5p
ggplot(data_5p, aes(x = Position, y = Entropy)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = c(6,12), linetype="dashed") +
  annotate("rect", xmin=6, xmax=12, ymin=-Inf, ymax=Inf, alpha=0.2) +
  theme_minimal() +
  labs(title = "Shannon Entropy per Position (miR-548-5p sequences)",
       x = "Position",
       y = "Entropy")
