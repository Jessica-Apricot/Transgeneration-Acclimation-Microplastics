#install.packages("BiocManager")
#BiocManager::install("limma")

library(BiocManager)
library(limma)
library(Rsubread)

# List all ZIP files in the current directory
zip_files <- list.files(pattern = "\\.zip$")

# Loop through and unzip each file
for (file in zip_files) {
  unzip(file)
}


# List all "summary.txt" files in subdirectories
summary_files <- list.files(path = ".", pattern = "summary.txt$", recursive = TRUE, full.names = TRUE)

# Read and concatenate all files
summary_content <- sapply(summary_files, readLines, simplify = FALSE)

# Flatten into a single vector
summary_content <- unlist(summary_content)

# Write to the output file
writeLines(summary_content, "~/Microplastics_RNA-seq/Microplastics_Raw_data/fastq/Brain/QC/fastqc_summaries.txt")


#######################################

# Define your directories for each tissue type
brain_dir <- "C:/Users/Jessi/Documents/Microplastics_RNA-seq/Microplastics_Raw_data/fastq/Brain/Hisat2"
gonad_dir <- "C:/Users/Jessi/Documents/Microplastics_RNA-seq/Microplastics_Raw_data/fastq/Gonad/Hisat2"
liver_dir <- "C:/Users/Jessi/Documents/Microplastics_RNA-seq/Microplastics_Raw_data/fastq/Liver/hisat2"


#BAM Files list
bam_files <- c(
  list.files(path = brain_dir, pattern = ".*\\.bam$", full.names = TRUE),
  list.files(path = gonad_dir, pattern = ".*\\.bam$", full.names = TRUE),
  list.files(path = liver_dir, pattern = ".*\\.bam$", full.names = TRUE)
)


bam_files_str <- paste(shQuote(bam_files), collapse = " ")
gtf_file <- "C:/Users/Jessi/Documents/Microplastics_RNA-seq/Genome/GCF_000002035.6_GRCz11_genomic.gtf"
gtf_file_quoted <- shQuote(gtf_file)

fc_results <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  nthreads = 2
)


#Make a table of feature counts
gene_info <- fc_results$annotation[, c("GeneID", "Chr", "Start", "End", "Strand", "Length")]

gene_info$Start <- sapply(strsplit(gene_info$Start, ";"), `[`, 1)
gene_info$End <- sapply(strsplit(gene_info$End, ";"), `[`, 1)

gene_info$Start <- as.numeric(gene_info$Start)
gene_info$End <- as.numeric(gene_info$End)

counts <- fc_results$counts
sample_names <- gsub(".bam$", "", basename(bam_files))
counts <- as.matrix(fc_results$counts)
colnames(counts) <- sample_names
count_data <- cbind(gene_info, counts)


write.table(gene_info, file = "gene_info.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(count_data, file = "C:/Users/Jessi/Documents/Microplastics_RNA-seq/Results/count_data.csv", row.names = FALSE)


#Renaming Samples
names(count_data)[7:78] <- c("MP_F0_M1_B", "D_F0_M1_B", "L_F0_M1_B", "H_F1_M3_B", "H_F1_M3_B.2.", "D_F0_M2_B", "D_F1_M3_B", "MPD_F0_M1_B", "D_F0_M3_B", "L_F0_M2_B", "L_F0_M3_B", "H_F0_M2_B", "H_F0_M1_B", "H_F0_M3_B", "MP_F0_M3_B", "MP_F0_M2_B", "MPD_F0_M2_B", "MPD_F0_M3_B", "D_F1_M1_B", "D_F1_M1_B.2.", "H_F1_M2_B", "MP_F0_M4_B", "H_F0_M4_B", "D_F0_M4_B", "D_F0_M2_G", "L_F0_M2_G", "L_F0_M3_G", "H_F0_M2_G","D_F1_M1_G", "D_F1_M1_G.2.", "D_F1_M3_G", "H_F1_M3_G", "H_F1_M3_G.2.", "H_F1_M2_G", "MP_F0_M1_G", "MPD_F0_M1_G", "H_F0_M1_G", "H_F0_M3_G", "MP_F0_M3_G", "MP_F0_M2_G", "MPD_F0_M2_G", "MPD_F0_M3_G", "D_F0_M3_G", "D_F0_M1_G", "L_F0_M1_G", "MP_F0_M4_G", "H_F0_M4_G", "D_F0_M4_G", "L_F0_M1_L", "L_F0_M2_L", "L_F0_M3_L", "H_F0_M1_L", "MP_F0_M2_L", "D_F1_M3_L", "H_F1_M3_L", "H_F1_M3_L.2.", "H_F1_M2_L", "MP_F0_M1_L", "MPD_F0_M1_L", "D_F0_M3_L", "D_F0_M2_L", "D_F0_M1_L", "H_F0_M2_L", "H_F0_M3_L", "MP_F0_M3_L", "MPD_F0_M2_L", "MPD_F0_M3_L", "D_F1_M1_L", "D_F1_M1_L.2.", "MP_F0_M4_L", "H_F0_M4_L", "D_F0_M4_L" )

All_Counts <- count_data[,7:78]

Brain_Counts <- count_data[,7:30]
Gonad_Counts <- count_data[,31:54]
Liver_Counts <- count_data[,55:78]

write.table(All_Counts, file = "All_Counts.txt", sep = "\t", quote = FALSE, col.names = NA)
write.table(Brain_Counts, file = "Brain_Counts.txt", sep = "\t", quote = FALSE, col.names = NA)
write.table(Gonad_Counts, file = "Gonad_Counts.txt", sep = "\t", quote = FALSE, col.names = NA)
write.table(Liver_Counts, file = "Liver_Counts.txt", sep = "\t", quote = FALSE, col.names = NA)

rm(summary_files, brain_dir, gonad_dir, liver_dir, bam_files, gtf_file, bam_files_str, gtf_file_quoted)



zero <- rownames(All_Counts)[rowSums(All_Counts == 0) == ncol(All_Counts)]
length(zero)
#WE have 2484 with zero expression across all samples. 
colSums(All_Counts)


colSums(All_Counts) %>% barplot(., las = 3, ylab = "Reads Mapped Per Sample")

boxplot(as.matrix(All_Counts) ~ col(All_Counts), 
        ylab = "Counts", xlab = "Samples", names = colnames(All_Counts))


#Brain Graphs
colSums(Brain_Counts)
colSums(Brain_Counts) %>% barplot(., las = 3, ylab = "Reads Mapped Per Sample", cex.names = 0.7, main = "Brain Samples")

#Plus 0.5 to avoid log(0) errors
logBrain = log2(as.matrix(Brain_Counts) + 0.5)
boxplot(as.matrix(logBrain) ~ col(Brain_Counts), las=3, cex.axis = 0.4,
        ylab = "Counts", xlab = "Samples", names = colnames(Brain_Counts), main = "Log counts of Brain samples")

#Gonad Graphs
colSums(Gonad_Counts)
colSums(Gonad_Counts) %>% barplot(., las = 3, ylab = "Reads Mapped Per Sample", cex.names = 0.7, main = "Gonad Samples")

logGonad = log2(as.matrix(Gonad_Counts) + 0.5)
boxplot(as.matrix(logGonad) ~ col(Gonad_Counts), las=3, cex.axis = 0.4,
        ylab = "Counts", xlab = "Samples", names = colnames(Gonad_Counts), main = "Log counts of Gonad samples")

#Liver Graphs
colSums(Liver_Counts)
colSums(Liver_Counts) %>% barplot(., las = 3, ylab = "Reads Mapped Per Sample", cex.names = 0.7, main = "Liver Samples")

logLiver = log2(as.matrix(Liver_Counts) + 0.5)
boxplot(as.matrix(logLiver) ~ col(Liver_Counts), las=3, cex.axis = 0.4,
        ylab = "Counts", xlab = "Samples", names = colnames(Liver_Counts), main = "Log counts of Liver samples")