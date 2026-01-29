#!/usr/bin/env Rscript
library("Rsamtools")
library("rtracklayer")
library("GenomicRanges")
library("devtools")
library("rstudioapi")
library("derfinder")
library("ggplot2")
library("caret")
library("reshape2")
library("data.table")
library("BiocParallel")



files <- list.files()                # Get all files in the directory
files <- files[endsWith(files, ".bam")]


cutoff <- 5
#all chromosomes of hg38 used, not using Y chromosomes 
#chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrM","chr1_KI270706v1_random","chr1_KI270707v1_random","chr1_KI270708v1_random","chr1_KI270709v1_random","chr1_KI270710v1_random","chr1_KI270711v1_random","chr1_KI270712v1_random","chr1_KI270713v1_random","chr1_KI270714v1_random","chr2_KI270715v1_random","chr2_KI270716v1_random","chr3_GL000221v1_random","chr4_GL000008v2_random","chr5_GL000208v1_random","chr9_KI270717v1_random","chr9_KI270718v1_random","chr9_KI270719v1_random","chr9_KI270720v1_random","chr11_KI270721v1_random","chr14_GL000009v2_random","chr14_GL000225v1_random","chr14_KI270722v1_random","chr14_GL000194v1_random","chr14_KI270723v1_random","chr14_KI270724v1_random","chr14_KI270725v1_random","chr14_KI270726v1_random","chr15_KI270727v1_random","chr16_KI270728v1_random","chr17_GL000205v2_random","chr17_KI270729v1_random","chr17_KI270730v1_random","chr22_KI270731v1_random","chr22_KI270732v1_random","chr22_KI270733v1_random","chr22_KI270734v1_random","chr22_KI270735v1_random","chr22_KI270736v1_random","chr22_KI270737v1_random","chr22_KI270738v1_random","chr22_KI270739v1_random","chrUn_KI270302v1","chrUn_KI270304v1","chrUn_KI270303v1","chrUn_KI270305v1","chrUn_KI270322v1","chrUn_KI270320v1","chrUn_KI270310v1","chrUn_KI270316v1","chrUn_KI270315v1","chrUn_KI270312v1","chrUn_KI270311v1","chrUn_KI270317v1","chrUn_KI270412v1","chrUn_KI270411v1","chrUn_KI270414v1","chrUn_KI270419v1","chrUn_KI270418v1","chrUn_KI270420v1","chrUn_KI270424v1","chrUn_KI270417v1","chrUn_KI270422v1","chrUn_KI270423v1","chrUn_KI270425v1","chrUn_KI270429v1","chrUn_KI270442v1","chrUn_KI270466v1","chrUn_KI270465v1","chrUn_KI270467v1","chrUn_KI270435v1","chrUn_KI270438v1","chrUn_KI270468v1","chrUn_KI270510v1","chrUn_KI270509v1","chrUn_KI270518v1","chrUn_KI270508v1","chrUn_KI270516v1","chrUn_KI270512v1","chrUn_KI270519v1","chrUn_KI270522v1","chrUn_KI270511v1","chrUn_KI270515v1","chrUn_KI270507v1","chrUn_KI270517v1","chrUn_KI270529v1","chrUn_KI270528v1","chrUn_KI270530v1","chrUn_KI270539v1","chrUn_KI270538v1","chrUn_KI270544v1","chrUn_KI270548v1","chrUn_KI270583v1","chrUn_KI270587v1","chrUn_KI270580v1","chrUn_KI270581v1","chrUn_KI270579v1","chrUn_KI270589v1","chrUn_KI270590v1","chrUn_KI270584v1","chrUn_KI270582v1","chrUn_KI270588v1","chrUn_KI270593v1","chrUn_KI270591v1","chrUn_KI270330v1","chrUn_KI270329v1","chrUn_KI270334v1","chrUn_KI270333v1","chrUn_KI270335v1","chrUn_KI270338v1","chrUn_KI270340v1","chrUn_KI270336v1","chrUn_KI270337v1","chrUn_KI270363v1","chrUn_KI270364v1","chrUn_KI270362v1","chrUn_KI270366v1","chrUn_KI270378v1","chrUn_KI270379v1","chrUn_KI270389v1","chrUn_KI270390v1","chrUn_KI270387v1","chrUn_KI270395v1","chrUn_KI270396v1","chrUn_KI270388v1","chrUn_KI270394v1","chrUn_KI270386v1","chrUn_KI270391v1","chrUn_KI270383v1","chrUn_KI270393v1","chrUn_KI270384v1","chrUn_KI270392v1","chrUn_KI270381v1","chrUn_KI270385v1","chrUn_KI270382v1","chrUn_KI270376v1","chrUn_KI270374v1","chrUn_KI270372v1","chrUn_KI270373v1","chrUn_KI270375v1","chrUn_KI270371v1","chrUn_KI270448v1","chrUn_KI270521v1","chrUn_GL000195v1","chrUn_GL000219v1","chrUn_GL000220v1","chrUn_GL000224v1","chrUn_KI270741v1","chrUn_GL000226v1","chrUn_GL000213v1","chrUn_KI270743v1","chrUn_KI270744v1","chrUn_KI270745v1","chrUn_KI270746v1","chrUn_KI270747v1","chrUn_KI270748v1","chrUn_KI270749v1","chrUn_KI270750v1","chrUn_KI270751v1","chrUn_KI270752v1","chrUn_KI270753v1","chrUn_KI270754v1","chrUn_KI270755v1","chrUn_KI270756v1","chrUn_KI270757v1","chrUn_GL000214v1","chrUn_KI270742v1","chrUn_GL000216v2","chrUn_GL000218v1")
chrs = readLines("chrName.txt")
Rlen = 22
#/New_data/Chromosome_calling/new_version/gencodeV33_pluspirnadb1_7_6_plusmirbase21_tRNAscan_MINT_hg38_V2.gff3
#gffdata<-import.gff("/New_data/Chromosome_calling/new_version/gencodeV33_pluspirnadb1_7_6_plusmirbase21_tRNAscan_MINT_hg38_V2.gff3")
gffdata<-import.gff("subset.gtf1")

#gffdata rna central
gffdata2<-import.gff("subset.gtf2")
names(mcols(gffdata2)) <- c("source", "type", "phase", "Name", "type.1", "databases", "ID", "source.1", "Parent", "providing_databases")

#gffdata <- within(gffdata, rm(score))
mcols(gffdata)\$score <- NULL
seqlevelsStyle(gffdata) <- "UCSC"

#gffdata2 <- within(gffdata2, rm(score))
mcols(gffdata2)\$score <- NULL
seqlevelsStyle(gffdata2) <- "UCSC"

# We recommend using a single feature-type annotation to simplify interpretation, in this case gene, optional step
gffdata<- gffdata[gffdata\$type == "gene",]


options(expressions = 500000)

Get_Annotated_Matrix<-function(chrs,cutoff, bams, Len, anno){
  
  #Create fullCoverage object for derfinder
  print("Creating Full Coverage Object...")
  fullCov <- fullCoverage(files = bams, chrs = chrs, verbose = T, mc.cores = 10, mapqFilter=255)
  #filteredCov <- lapply(fullCov, filterData, cutoff = cutoff)
  #save.image(file = "/mnt/acri7/tmp.RData")
  #rm(fullCov)
  #print('dimentions of coverage object')
  #print(dim(filteredCov))
  #Get library sizes
  TM<-vector(mode="integer", length = 0)
  
  for (i in seq(1,length(bams))){
    Tmapped<-getTotalMapped(bams[[i]],chrs = chrs)
    TM<-c(TM,Tmapped[1])
  }
  
  #Get Expressed-region level counts for fullCoverage object
  print("Extracting Count matrix from derfinder object...")
  regionMat <- list()
  regionMat <- regionMatrix(fullCov, cutoff = cutoff, L = Len, verbose = FALSE, targetSize = mean(TM), totalMapped = TM)
  #save(regionMat, file = "/mnt/acri7/tmp_regionmat.Rda")
  #Extract data from RegionMatrix Object
  GRL<-GRangesList()
  
  #Create GrangesList object containing the information for each contig. 
  for (i in chrs) {
    values(regionMat[[i]]\$regions)<-regionMat[[i]]\$coverageMatrix
    GRL<-c(GRL,GRangesList(regionMat[[i]]\$regions))
  }
  
  #Unlist this Grange list to get an unlisted Granges object. 
  UL<-unlist(GRL)

  #find overlaps between the expressed regions and the target annotation
  hits1<-findOverlaps(UL,gffdata)
  hits1_1 <-findOverlaps(UL,gffdata2)
  #custom function to replace a vector of redundant strings with only unique values
  pasteu<-function(x){
    paste(unique(unlist(strsplit(x,";"))),collapse =';')
  }
  
  #convert to dataframe
  print("Finding Overlaps between annotation and Expressed Regions...")
  hits1<-data.frame(hits1)
  hits1<-as.data.table(hits1)
  hits1_1 <- data.frame(hits1_1)
  hits1_1 <- as.data.table(hits1_1)
  
  #isolate gene IDs, RNA types and % coverage from the GTF annotation with associated indices
  hits1\$id<-NA
  hits1\$type<-NA
  hits1\$perc_cov<-NA
  hits1\$id<-sapply(gffdata\$gene_name[(hits1\$subjectHits)], `[[`, 1)
  hits1\$type<-unlist(gffdata\$gene_type[(hits1\$subjectHits)])
  hits1\$perc_cov<-round(width(UL)[hits1\$queryHits]/width(gffdata)[(hits1\$subjectHits)],2)
  hits1_1\$id_rnaCentral <- NA
  hits1_1\$type_rnaCentral <- NA
  hits1_1\$perc_cov_rnaCentral <- NA
  hits1_1\$id_rnaCentral<-sapply(gffdata2\$Name[(hits1_1\$subjectHits)], `[[`, 1)
  hits1_1\$type_rnaCentral<-unlist(gffdata2\$type.1[(hits1_1\$subjectHits)])
  hits1_1\$perc_cov_rnaCentral<-round(width(UL)[hits1_1\$queryHits]/width(gffdata2)[(hits1_1\$subjectHits)],2)
  
  #Aggregate findOverlaps hits by query indices while reducing the overlapping annotations to get only unique IDs
  hits2<-aggregate(hits1\$id,list(hits1\$queryHits),paste,collapse=';')
  hits2_1<-aggregate(hits1_1\$id_rnaCentral,list(hits1_1\$queryHits),paste,collapse=';')
  reduced_ids<-lapply(hits2\$x,pasteu)
  reduced_ids<-unlist(reduced_ids)
  hits2\$x<-reduced_ids
  reduced_ids <- lapply(hits2_1\$x,pasteu)
  reduced_ids <- unlist(reduced_ids)
  hits2_1\$x <- reduced_ids
  
  hits3<-aggregate(hits1\$type,list(hits1\$queryHits),paste,collapse=';')
  reduced_types<-lapply(hits3\$x,pasteu)
  reduced_types<-unlist(reduced_types)
  hits3\$x<-reduced_types
  hits3_1 <- aggregate(hits1_1\$type_rnaCentral,list(hits1_1\$queryHits),paste,collapse=';')
  reduced_types <- lapply (hits3_1\$x, pasteu)
  reduced_types <- unlist(reduced_types)
  hits3_1\$x <- reduced_types
  
  
  hits4<-aggregate(hits1\$perc_cov,list(hits1\$queryHits),paste,collapse=';')
  reduced_cov<-lapply(hits4\$x,pasteu)
  reduced_cov<-unlist(reduced_cov)
  hits4\$x<-reduced_cov

  hits4_1<-aggregate(hits1_1\$perc_cov_rnaCentral,list(hits1_1\$queryHits),paste,collapse=';')
  reduced_cov<-lapply(hits4_1\$x,pasteu)
  reduced_cov<-unlist(reduced_cov)
  hits4_1\$x<-reduced_cov
  
  #Create a dataframe object from the annotated subset
  print("Building Final Tables")
  mt_dfr<-cbind(as.character(seqnames(UL)),start(UL)-1,end(UL))
  colnames(mt_dfr)<-c("contig","start","end")
  mt_dfr<-cbind(mt_dfr,values(UL))

  mt_dfr\$contig<-as.character(mt_dfr\$contig)
  mt_dfr\$start<-as.character(mt_dfr\$start)
  mt_dfr\$end<-as.character(mt_dfr\$end)

  #Append gene ids and RNA type info to expressed regions
  mt_dfr\$ids<-NA
  mt_dfr\$type<-NA
  mt_dfr\$perc_cov<-NA
  mt_dfr\$ids_rnaCentral <- NA
  mt_dfr\$type_rnaCentral <- NA
  mt_dfr\$perc_cov_rnaCentral <- NA
  mt_dfr\$ids[hits2\$Group.1]<-hits2\$x
  mt_dfr\$type[hits3\$Group.1]<-hits3\$x
  mt_dfr\$perc_cov[hits4\$Group.1]<-hits4\$x
  mt_dfr\$ids_rnaCentral[hits2_1\$Group.1] <- hits2_1\$x
  mt_dfr\$type_rnaCentral[hits3_1\$Group.1] <- hits3_1\$x
  mt_dfr\$perc_cov_rnaCentral[hits4_1\$Group.1] <- hits4_1\$x
  #Return Complete Object
  return(mt_dfr)
}

#Example of annotated Count Matrix with Gencode GFF

print(files)


df_1<- Get_Annotated_Matrix(chrs,cutoff,files,Rlen, gffdata)
save(df_1, file = "derfinder_out_parallel_unique.Rda")






# Get R version
r.version <- strsplit(version[['version.string']], ' ')[[1]][3]

# Get package versions
packages <- c(
	        "Rsamtools", "rtracklayer", "GenomicRanges", "devtools", "rstudioapi",
		  "derfinder", "ggplot2", "caret", "reshape2", "data.table", "BiocParallel"
		)

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
