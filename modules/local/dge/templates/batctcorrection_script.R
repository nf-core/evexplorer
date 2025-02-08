#!/usr/bin/env Rscript
library("sva")


sample_sheet <- read.csv("samplesheet.csv")
load("derfinder_out_parallel_unique.Rda")  # Replace with the actual file name

ls()  # Check the loaded objects
head(df_1)  # Inspect the dataframe


selected_samples <- samplesheet$sample


colnames(df_1) <- gsub("\\.bam$", "", colnames(df_1))
colnames(df_1) <- gsub("\\.", "-", colnames(df_1))


extracted_df_1 <- df_1[, colnames(df_1) %in% selected_samples]  # Subset the data

count_matrix <- ComBat_seq(as.matrix(extracted_df_1), batch=samplesheet$batch, group=NULL)

write.csv(count_matrix, file="count_matrix.csv", row.names=TRUE,quote=FALSE)

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
