library("SummarizedExperiment") # bioconductor package
library("DESeq2") # bioconductor package
library("IHW") # bioconductor package
library("biomaRt") # bioconductor package
library("apeglm") # bioconductor package
library("pheatmap") # CRAN package
library("RColorBrewer") # CRAN package
library("PCAtools") # bioconductor package
library(reshape2) # CRAN package
library(dplyr)
library(purrr)

setwd("/Users/shivgurjar/Documents/Terminal/BUCCALMUCOSA/counts")

# Get a list of all feature count files in the directory
file_paths <- list.files(pattern = "_featurecounts.txt", full.names = TRUE)

# Read data and create a list of data frames
countDataList <- map(file_paths, ~{
  # Read data with check.names set to FALSE
  data <- read.table(.x, header=TRUE, row.names=1, check.names = FALSE)
  data <- subset(data, select = -c(Chr, Start, End, Strand, Length))
  # Add a column for sample type (if needed)
  # data$SampleType <- "YourSampleType"
  return(data)
})

# Remove NULL entries (files that don't exist)
countDataList <- compact(countDataList)
print(countDataList)

# Combine data frames into a single data frame
combined_data <- do.call(cbind, countDataList)
print(combined_data)

# Remove ".bam" from column names
colnames(combined_data) <- gsub("\\.bam|/", "", colnames(combined_data))
colnames(combined_data) <- gsub("bam_", "", colnames(combined_data))

# Print the combined data
print(combined_data)
str(countDataList)


write.table(colData, file = "colData.txt", sep = "\t", quote = FALSE, row.names = TRUE)
colData <- data.frame(
  condition = rep("Control", 3),
  age = c(48, 50, 49.68),
  subsite = rep("Mouth Floor", 3),
  sex = c("male", "female", "male"),
  SRR_ID = c("SRR7415602", "SRR7415601", "SRR7415599")
)
# Additional metadata
additional_metadata <- data.frame(
  condition = rep("Tumor", 4),
  age = c(56, 53, 61, 66),
  subsite = rep("Mouth Floor", 4),
  sex = c("male", "male", "male", "male"),
  SRR_ID = c("SRR21141186 ", "SRR21141185 ", "SRR21141184 ", "SRR21141183 ")
)
# Combine colData and additional_metadata
colData <- rbind(colData, additional_metadata)
print(colData)


# Print dimensions of combined_data and colData
cat("Dimensions of combined_data:", dim(combined_data), "\n")
cat("Dimensions of colData:", dim(colData), "\n")
# ... rest of your code
# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = combined_data,
  colData = colData,
  design = ~ condition
)

#pre-filtering to keep only rows that have at least 1 reads total
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

# Run the DESeq
dds <- DESeq(dds)


res <- results(dds)
head(results(dds, tidy=TRUE))

summary(res) #summary of results

# variance stabilize data
rld <- rlog(dds)
head(assay(rld))

# Exploratory data analysis
# Plot top 500 most variable genes
vars <- rowVars(assay(rld))
vars <- sort(vars, decreasing=TRUE)
vars <- vars[1:500]
topVars <- assay(rld)[names(vars),]



# Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
tiff("MAplot.tiff")

# Set Cook's distance cutoff
res_cook <- results(dds, cooksCutoff=0.01)

plotMA(res_cook, ylim=c(-7,7))
abline(h=c(-1,1), col="blue")
dev.off()
sum(res$padj < 0.1, na.rm=TRUE)

#summary(res05)

#sum(res05$padj < 0.1, na.rm=TRUE)

# Use independent filtering with a specified threshold
res_filtered <- results(dds, independentFiltering=TRUE, alpha=0.1)

# Display the summary for the Cook's distance cutoff
summary(res_cook)

# Display the summary for the filtered results
summary(res_filtered)
