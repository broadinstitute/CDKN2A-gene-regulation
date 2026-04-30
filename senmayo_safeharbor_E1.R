## ----------------------------
## Libraries
## ----------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

## ----------------------------
## Load Seurat object
## ----------------------------
screen.rna <- readRDS("data/screen.rna.rds")

## ----------------------------
## Guide columns
## ----------------------------
guide_cols <- colnames(screen.rna@meta.data)[7:730]
safe_harbor_guides <- guide_cols[grepl("^safe\\.harbor", guide_cols)]

## Target guides (your high MOI gene of interest)
guide_names <- colnames(screen.rna@meta.data)[79:82]
threshold <- 10

## ----------------------------
## Filter cells by num_guides <= 50
## ----------------------------
binary_gRNA_matrix <- ifelse(screen.rna@meta.data[, guide_cols] >= threshold, 1, 0)
num_guides <- rowSums(binary_gRNA_matrix)
cells_to_keep <- rownames(screen.rna@meta.data)[num_guides <= 50]

print(paste("Cells before num_guides filter:", ncol(screen.rna)))
print(paste("Cells after num_guides <= 50 filter:", length(cells_to_keep)))

# Subset Seurat object
screen.rna <- subset(screen.rna, cells = cells_to_keep)

## ----------------------------
## Calculate SenMayo signature
## ----------------------------
senmayo_genes <- read.table(
  "SenMayo_signature.txt",
  stringsAsFactors = FALSE
)[,1]

# Keep only genes present in the dataset
senmayo_genes <- intersect(senmayo_genes, rownames(screen.rna))
print(paste("SenMayo genes in dataset:", length(senmayo_genes)))

# Calculate mean expression of SenMayo genes per cell
senmayo_expr <- colMeans(as.matrix(screen.rna@assays$RNA@data[senmayo_genes, ]))
screen.rna$SenMayo_score <- senmayo_expr

## ----------------------------
## Identify cells with target vs safe harbor guides
## ----------------------------
# Cells with your target guides
cells_with_target <- rownames(screen.rna@meta.data)[
  rowSums(screen.rna@meta.data[, guide_names] >= threshold) > 0
]

# Cells with safe harbor guides
cells_with_safeharbor <- rownames(screen.rna@meta.data)[
  rowSums(screen.rna@meta.data[, safe_harbor_guides] >= threshold) > 0
]

# Cells with BOTH target and safe harbor - EXCLUDE THESE
cells_with_both <- intersect(cells_with_target, cells_with_safeharbor)

print(paste("Cells with target guides:", length(cells_with_target)))
print(paste("Cells with safe harbor guides:", length(cells_with_safeharbor)))
print(paste("Cells with BOTH (to be excluded):", length(cells_with_both)))

# Remove cells with both
cells_with_target_only <- setdiff(cells_with_target, cells_with_both)
cells_with_safeharbor_only <- setdiff(cells_with_safeharbor, cells_with_both)

print(paste("Cells with target ONLY:", length(cells_with_target_only)))
print(paste("Cells with safe harbor ONLY:", length(cells_with_safeharbor_only)))

# Create comparison column
screen.rna$comparison_group <- ifelse(
  rownames(screen.rna@meta.data) %in% cells_with_target_only, "Target",
  ifelse(rownames(screen.rna@meta.data) %in% cells_with_safeharbor_only, "SafeHarbor", "Other")
)

## ----------------------------
## Calculate num_guides for each cell
## ----------------------------
screen.rna$num_guides <- rowSums(screen.rna@meta.data[, guide_cols] >= threshold)

## ----------------------------
## Subset for comparison
## ----------------------------
df_compare <- screen.rna@meta.data[
  screen.rna@meta.data$comparison_group %in% c("SafeHarbor", "Target"),
]

print(paste("Cells in final comparison:", nrow(df_compare)))
print(table(df_compare$comparison_group))

## ----------------------------
## Linear model controlling for num_guides
## ----------------------------
fit <- lm(SenMayo_score ~ comparison_group + num_guides, data = df_compare)
summary(fit)

## Extract coefficient and p-value for Target effect
coef_summary <- summary(fit)$coefficients
print(coef_summary)

## ----------------------------
## Calculate residuals
## ----------------------------
df_compare$SenMayo_resid <- resid(fit)

## ----------------------------
## Plot 1: Residualized SenMayo score
## ----------------------------
p1 <- ggplot(df_compare,
             aes(x = comparison_group,
                 y = SenMayo_resid,
                 fill = comparison_group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  labs(
    title = "SenMayo Score (num_guides regressed out)",
    subtitle = paste0("Target n=", sum(df_compare$comparison_group == "Target"),
                      ", SafeHarbor n=", sum(df_compare$comparison_group == "SafeHarbor")),
    x = NULL,
    y = "SenMayo residual"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave(
  filename = "SenMayo_residual_violin.png",
  plot = p1,
  width = 5,
  height = 4,
  dpi = 300
)

## ----------------------------
## Plot 2: Sanity check - num_guides effect
## ----------------------------
p2 <- ggplot(df_compare,
             aes(x = num_guides,
                 y = SenMayo_score,
                 color = comparison_group)) +
  geom_point(alpha = 0.3, size = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "SenMayo Score vs Number of Guides",
    x = "Number of guides",
    y = "SenMayo score"
  ) +
  theme_classic(base_size = 14)

ggsave(
  filename = "SenMayo_vs_num_guides.png",
  plot = p2,
  width = 6,
  height = 4,
  dpi = 300
)

## ----------------------------
## Plot 3: Distribution of num_guides by group
## ----------------------------
p3 <- ggplot(df_compare,
             aes(x = comparison_group,
                 y = num_guides,
                 fill = comparison_group)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  labs(
    title = "Number of Guides by Group",
    x = NULL,
    y = "Number of guides"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave(
  filename = "num_guides_distribution.png",
  plot = p3,
  width = 5,
  height = 4,
  dpi = 300
)

## ----------------------------
## Statistical tests
## ----------------------------
# Wilcoxon test on residuals (non-parametric, doesn't assume normality)
wilcox_result <- wilcox.test(SenMayo_resid ~ comparison_group, data = df_compare)
print("Wilcoxon test on residuals:")
print(wilcox_result)

# T-test on residuals (parametric alternative)
t_result <- t.test(SenMayo_resid ~ comparison_group, data = df_compare)
print("T-test on residuals:")
print(t_result)

# Check if num_guides differs between groups (sanity check)
numguides_test <- wilcox.test(num_guides ~ comparison_group, data = df_compare)
print("Wilcoxon test - num_guides difference between groups:")
print(numguides_test)
