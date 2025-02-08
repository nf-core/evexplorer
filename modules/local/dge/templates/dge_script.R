library(edgeR)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)


# Read count matrix
count_matrix <- read.csv("count_matrix.csv", row.names = 1, header = TRUE)

# Read sample metadata
samplesheet <- read.csv("comboseq_samplesheet.csv", stringsAsFactors = FALSE)

# Extract relevant columns (batch & condition)
design <- samplesheet[, c("sample", "batch", "condition")]
row.names(design) <- design$sample



# changing . in column names with hyphen
colnames(count_matrix) <- gsub("\\.", "-", colnames(count_matrix))

# Ensure column names of count matrix match sample names in design
count_matrix <- count_matrix[, design$sample]



# Convert batch and condition to factors
design$batch <- as.factor(design$batch)
design$condition <- as.factor(design$condition)



# Create a DGEList object
dge <- DGEList(counts = count_matrix, group = design$condition)

# Filter out lowly expressed genes (genes with very low counts)
keep <- filterByExpr(dge, group = design$condition)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalize using TMM method
dge <- calcNormFactors(dge)


# Create model matrix
design_matrix <- model.matrix(~ batch + condition, data = design)

# Estimate dispersion
dge <- estimateDisp(dge, design_matrix)

# Fit the model
fit <- glmFit(dge, design_matrix)

# Perform differential expression analysis
lrt <- glmLRT(fit, coef = "conditionnormal")  # Change as per your levels

# Extract top DE genes
top_genes <- topTags(lrt, n=50)$table

head(top_genes)

pdf("DE_analysis_plots.pdf", width = 10, height = 7)

# 5.1. MDS Plot (Multidimensional Scaling)
plotMDS(dge, col = as.numeric(design$condition), main = "MDS Plot")

# 5.2. PCA Plot
pca_res <- prcomp(t(log2(dge$counts + 1)))
pca_data <- as.data.frame(pca_res$x)
pca_data$condition <- design$condition

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 4) +
    labs(title = "PCA Plot") +
    theme_minimal()
print(ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 4) +
    labs(title = "PCA Plot") +
    theme_minimal())

# 5.3. MA Plot
plotMD(lrt, main = "MA Plot", ylim = c(-5, 5))

# 5.4. Volcano Plot
volcano_data <- data.frame(logFC = lrt$table$logFC, 
                           pval = -log10(lrt$table$PValue))
ggplot(volcano_data, aes(x = logFC, y = pval)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(title = "Volcano Plot", x = "Log Fold Change", y = "-log10 P-value") +
    theme_minimal()
print(ggplot(volcano_data, aes(x = logFC, y = pval)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(title = "Volcano Plot", x = "Log Fold Change", y = "-log10 P-value") +
    theme_minimal())

# 5.5. Heatmap of Top DE Genes
top_gene_ids <- rownames(top_genes)
pheatmap(log2(dge$counts[top_gene_ids, ] + 1), 
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         main = "Heatmap of Top DE Genes")

dev.off()

