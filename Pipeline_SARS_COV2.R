# Download reference genome from human genome 38 from bowtie2 

# Download annotation file from the database 
hsa.gff3 file 

################ SRA FILES TO FASTQC FILES ###############

## Path to sratoolkit module sra-tool/2.9.6
sra_path = "/N/soft/rhel7/sra-toolkit/sra-tool/2.9.6/bin/srapath

## Directory where the sra files are stored 
data_dir = "~/NGS/Assignment2/"
## Selecting all the files with .fastq format 
files = list.files(data_dir, pattern = ".sra$", full.names = TRUE)

## Output directory for aligned reads 
sra_fastq <- "~/NGS/Assignment2/SRROfastq" 


## Iterating through all the files and converting .sra to .fastq using fastq-dump
for (file in files) {
  output_file = file.path(aligned_read, basename(tools::file_path_sans_ext(file)))
  system(paste("fastq-dump", sra_file, "--outdir", fastq_output))
}

############# QUALITY CONTROL ##############
version 0.11.9
fastqc_path = "/N/soft/rhel7/fastqc/0.11.9/fastqc

for (file in fastq_files) {
  fastqc_cmd <- paste("fastqc", file)
  systemfastqc_path 
}


################ TRIMMING #####################

## Path to trimmoatic module 0.36
trim_path = "/N/soft/rhel7/trimmomatic/0.36/trimmomatic    

## Directory where the fastq files are stored 
data_dir = "~/NGS/Assignment2/SRROfastq"
## Selecting all the files with .fastq format 
files = list.files(data_dir, pattern = ".fastq$", full.names = TRUE)

## Output directory for aligned reads 
trim_files = "~/NGS/Assignment2/SRROfastq/Trim"

## Iterating through all the files 
for (file in files) {
  output_file = file.path(trim_files, basename(tools::file_path_sans_ext(file)))
  system(paste(trim_path, "SE", files, trim_files, "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10;"))
}



#################### ALIGNMENT ##################

## Path to the human genome reference downloaded from bowtie website - data is in .ht2 format 
genome_index = "~/NGS/Assignment2/SRROfastq/GRCh38_noalt_as"

## Path to bowtie2 module
bowtie2_path = "/N/soft/rhel7/bowtie2/2.4.2/bowtie2"

## Directory where the fastq files are stored 
data_dir = "~/NGS/Assignment2/SRROfastq/Trim"
## Selecting all the files with .fastq format 
files = list.files(data_dir, pattern = ".fastq$", full.names = TRUE)

## Output directory for aligned reads 
aligned_read = "~/NGS/Assignment2/SRROfastq/Alignment_output"

## Iterating through all the files 
for (file in files) {
  output_file = file.path(aligned_read, basename(tools::file_path_sans_ext(file)))
  system(paste(bowtie2_path, "-x", genome_index, "-U", files, "-S", output_file))
}



################## CONVERTING SAM TO BAM ###################
module load samtools/1.9

samtools view -bS 72_Trim-o 72.bam

# Sorting 
samtools view -bS 72.bam > 72.bam 


############### Feature counts #############
module load subread /2.0.6

featureCounts -t miRNA,miRNA_primary_transcript -g "Name" -a hsa.gff3 -o countsX.txt trim_X.bam


##### Reading Counts in R #########
SRR22269883 = read.table("counts_83.txt")
SRR22269883 = SRR22269883[,c(1,7)]
colnames(SRR22269883) = SRR22269883[1,]
SRR22269883 = SRR22269883[-1,]


SRR22269882 = read.table("counts_82.txt")
SRR22269882 = SRR22269882[,c(1,7)]
colnames(SRR22269882) = SRR22269882[1,]
SRR22269882 = SRR22269882[-1,]

.....

### Merging all the files 

counts_pre = merge(merge(merge(merge(SRR22269872, SRR22269873),SRR22269874),SRR22269875),SRR22269876)
counts_pre = merge(merge(merge(merge(merge(merge(counts_pre, SRR22269877),SRR22269878),SRR22269879),SRR22269880),SRR22269881), SRR22269882)

rownames(counts_pre) = counts_pre$Geneid

### Metadata is downloaded from SRA website ####

# Sample data frames (replace these with your actual data frames)
SARScov1_24 <- data.frame(SARScov1_24)
SARScov2_24 <- data.frame(SARScov2_24)
SARScov3_24 <- data.frame(SARScov3_24)

counts_pre <- merge(merge(merge(merge(SRR22269872,SRR22269873),SRR22269874),SRR22269875),SRR22269876)

counts_pre <- merge(merge(merge(merge(merge(counts_pre,SRR22269877),SRR22269878),SRR22269879),SRR22269880),SRR22269881)
counts_pre <- merge(merge(counts_pre,SRR22269882),SRR22269883)

rownames(counts_pre) <- counts_pre$Geneid
counts_pre$Geneid <- NULL

colnames(counts_pre) <- c('SRR22269872','SRR22269873','SRR22269874','SRR22269875','SRR22269876','SRR22269877','SRR22269878','SRR22269879','SRR22269880',
                          'SRR22269881','SRR22269882','SRR22269883')


counts_24_72 <- counts_pre[,which(colnames(counts_pre) %in% rownames(metadata_24_72))]

# List of data frames to merge
data_frames_to_merge <- list(SARScov1_24, SARScov2_24, SARScov3_24)

# Merge by the 'Geneid' column using Reduce
condition_24 <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), data_frames_to_merge)

# Print the merged data frame
print(condition_24)

# Sample data frames (replace these with your actual data frames)
SARScov1_72 <- data.frame(SARScov1_72)
SARScov2_72 <- data.frame(SARScov2_72)
SARScov3_72 <- data.frame(SARScov3_72)



# List of data frames to merge
data_frames_to_merge <- list(SARScov1_72, SARScov2_72, SARScov3_72)

# Merge by the 'Geneid' column using Reduce
condition_72 <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), data_frames_to_merge)

# Print the merged data frame
print(condition_72)

condition_72$condition <- c(rep("SARS-CoV-2 ", nrow(condition_72)))
condition_24$condition <- c(rep("SARS-CoV-2 ", nrow(condition_24)))




# Add a column for time
condition_72$time <- c(rep("72H", nrow(condition_72)))
condition_24$time <- c(rep("24H", nrow(condition_24)))

merged_counts <- merge(condition_24, condition_72, by = "Geneid")

# Assuming countMatrix is your count data matrix and merged_counts is your merged data frame
# Assuming countMatrix is your count data matrix and merged_counts is your merged data frame

dds = DESeqDataSetFromMatrix(
  countData = merged_counts,
  colData = DataFrame(Geneid = merged_counts$Geneid, time = merged_counts$time.x)
)

print(dim(merged_counts))

DESeq(
  "input" = merged_counts,
  "colCondition" = condition_24,
  "colTime" = time,
  "design" = ""
)


colData = merged_counts[, c("Geneid", "time.x")]
countData = countMatrix[, -1]
tcountData <- t(countData)

dim(countData)
dim(colData)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = countData[, -1],  # Exclude the first column (Geneid)
  colData = countData[, c("Geneid", "time.x", "time.y")],
  design = ~ time.x + time.y
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Contrast: Condition at 72H vs 24H
results_72vs24 <- results(dds, contrast = c("time", "72H", "24H"))

# Extract significant DEGs
significant_genes_72vs24 <- subset(results_72vs24, padj < 0.05)

# Print the significant DEGs
print(significant_genes_72vs24)






