#!/usr/bin/env Rscript

library("ggplot2")
library("pheatmap")
library("ggrepel")
library("RColorBrewer")
library("edgeR")

# Read count matrix
count_matrix <- read.csv("count_matrix.csv", row.names = 1, header = TRUE,check.names=FALSE)

clean_condition <- gsub("\\\\[|\\\\]", "", '$condition')


parts_condition <- unlist(strsplit(clean_condition, ",\\\\s*"))

# Turn into a matrix with 2 columns
mat <- matrix(parts_condition, ncol = 3, byrow = TRUE)


# Convert to dataframe
condition_info <- as.data.frame(mat, stringsAsFactors = FALSE)


# Set rownames from first column and drop that column
rownames(condition_info) <- condition_info[[1]]
design <- condition_info[ , -1, drop = FALSE]

colnames(design) <- c("batch","condition")


# Ensure column names of count matrix match sample names in design
count_matrix <- count_matrix[, rownames(design)]


# Convert batch and condition to factors
design\$batch <- as.factor(design\$batch)
design\$condition <- relevel(as.factor(design\$condition), ref = "normal")



# Ensure sample order matches
count_matrix <- count_matrix[, rownames(design)]

# Create a DGEList object
dge <- DGEList(counts = count_matrix, group = design\$condition)

# Filter out lowly expressed genes (genes with very low counts)
keep <- filterByExpr(dge, group = design\$condition)
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
lrt <- glmLRT(fit, coef = "conditiondisease")  # Change as per your levels

# Extract top DE genes
top_genes <- topTags(lrt, n=50)\$table

head(top_genes)

# writing the results to a file
write.csv(topTags(lrt, n = Inf)\$table,"differential_rna_expression.csv",quote=F)

pdf("DE_analysis_plots.pdf", width = 10, height = 7)


# 5.1. MDS Plot (Multidimensional Scaling)
plotMDS(dge, col = as.numeric(design\$condition), main = "MDS Plot")


# CPM for PCA
logCPM <- cpm(dge, log=TRUE, prior.count=1)

# 5.2. PCA Plot
pca_res <- prcomp(t(logCPM))
pca_data <- as.data.frame(pca_res\$x)
pca_data\$condition <- design\$condition

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
volcano_data <- data.frame(logFC = lrt\$table\$logFC, 
                           pval = -log10(lrt\$table\$PValue))
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
pheatmap(logCPM[top_gene_ids, ], 
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         main = "Heatmap of Top DE Genes")

dev.off()


# Get R version
r.version <- strsplit(version[['version.string']], ' ')[[1]][3]

# Get package versions
packages <- c("ggplot2","pheatmap","ggrepel","RColorBrewer","edgeR")

package.versions <- sapply(packages, function(pkg) {
				                                        tryCatch(as.character(packageVersion(pkg)), error = function(e) "NA")
									                })


# Write versions to YAML
writeLines(
	                c(
			                     '"${task.process}":',
					                            paste('    r-base:', r.version),
					                            sapply(names(package.versions), function(pkg) {
										                                               paste0('    r-', tolower(pkg), ': ', package.versions[[pkg]])
															                                                       })
								                             ),
	                'versions.yml'
			             )
