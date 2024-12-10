
# Load ggplot2 for plotting
library(ggplot2)

# Load data for ≤100kb
data_le <- read.table("data/processed/metric_max_100kb.txt", col.names = c("Length", "GC"))
# Load data for >100kb
data_gt <- read.table("data/processed/metric_min_100kb.txt", col.names = c("Length", "GC"))

# Plot length distribution for ≤100kb
ggplot(data_le, aes(x = Length)) +
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  scale_x_log10() +
  labs(title = "Sequence Length Distribution (≤ 100kb)", x = "Log10(Length)", y = "Frequency")

ggsave("data/output/length_histogram_max_100kb.png")

# Plot length distribution for >100kb
ggplot(data_gt, aes(x = Length)) +
  geom_histogram(bins = 30, fill = "red", color = "black") +
  scale_x_log10() +
  labs(title = "Sequence Length Distribution (> 100kb)", x = "Log10(Length)", y = "Frequency")

ggsave("data/output/length_histogram_min_100kb.png")



# Plot GC% distribution for ≤100kb
ggplot(data_le, aes(x = GC)) +
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  labs(title = "GC% Distribution (≤ 100kb)", x = "GC%", y = "Frequency")

ggsave("data/output/gc_histogram_max_100kb.png")


# Plot GC% distribution for >100kb
ggplot(data_gt, aes(x = GC)) +
  geom_histogram(bins = 30, fill = "red", color = "black") +
  labs(title = "GC% Distribution (> 100kb)", x = "GC%", y = "Frequency")

ggsave("data/output/gc_histogram_min_100kb.png")