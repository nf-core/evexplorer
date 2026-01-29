#!/usr/bin/env Rscript
library("sva")


load("derfinder_out_parallel_unique.Rda")


# Remove columns 2 and 3
colnames(df_1) <- gsub("\\\\.bam\$", "", colnames(df_1))
colnames(df_1) <- gsub("\\\\.", "-", colnames(df_1))


# Create the row names with the first colon between contig and start
rownames(df_1) <- paste(df_1\$contig, df_1\$start, df_1\$end, sep = ":")

# Replace the second colon with a hyphen using regular expression
rownames(df_1) <- gsub(":(\\\\d+):(\\\\d+)\$", ":\\\\1-\\\\2", rownames(df_1))



# Remove brackets from batch and condition
clean_batch <- gsub("\\\\[|\\\\]", "", '$batch')

clean_cond <- gsub("\\\\[|\\\\]", "", '$condition')





parts_batch <- unlist(strsplit(clean_batch, ",\\\\s*"))


parts_cond <- unlist(strsplit(clean_cond, ",\\\\s*"))



parts_condition <- parts_cond[!grepl("^[0-9]+\$", parts_cond)]



# Turn condition into a matrix with 2 columns
mat_cond <- matrix(parts_condition, ncol = 2, byrow = TRUE)


# Turn batch into a matrix with 2 columns
mat_batch <- matrix(parts_batch, ncol = 2, byrow = TRUE)




# Convert to dataframe
cond_info <- as.data.frame(mat_cond, stringsAsFactors = FALSE)

# Convert to dataframe
batch_info <- as.data.frame(mat_batch, stringsAsFactors = FALSE)




# Set rownames from first column and drop that column
rownames(batch_info) <- batch_info[[1]]
batch_info <- batch_info[ , -1, drop = FALSE]

colnames(batch_info) <- c("batch")


# Set rownames from first column and drop that column
rownames(cond_info) <- cond_info[[1]]
cond_info <- cond_info[ , -1, drop = FALSE]

colnames(cond_info) <- c("condition")

print(cond_info)

# filtering count matrix based on the sample names
extracted_df_1 <- df_1[, colnames(df_1) %in% rownames(batch_info)]  # Subset the data

count_matrix <- ComBat_seq(as.matrix(extracted_df_1), batch=batch_info\$batch, group=cond_info\$condition)

write.csv(count_matrix, file="count_matrix.csv",quote=FALSE)

# Get R version
r.version <- strsplit(version[['version.string']], ' ')[[1]][3]

# Get package versions
packages <- c( "sva")

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
