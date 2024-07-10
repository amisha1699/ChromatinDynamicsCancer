# Install and load the necessary packages
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")
BiocManager::install("biomaRt")
library(biomaRt)

library(org.Hs.eg.db)
library(AnnotationDbi)

# Assuming your DESeq results are stored in a variable named 'deseq_result'
# Extract the Ensemble gene IDs from your DESeq result
data <- read.csv("C:/Users/navee/OneDrive/Desktop/All-csv-files/combined-csv-files/deseq_result.csv")
head(data)
setwd("C:\\Users\\navee\\OneDrive\\Desktop")
getwd()
data_down <- read.delim("upgene-cutoff+1.txt",sep = "\n")
head(data_down)
# Assuming your DESeq results are stored in a variable named 'deseq_result'
# Extract the Ensemble gene IDs from your DESeq results
ensembl_ids <- data_down
ensembl_ids
head(ensembl_ids)
valid_keys <- keys(org.Hs.eg.db, keytype = "ENSEMBL")

# Check if any of the keys in ensembl_ids are not valid
invalid_keys <- ensembl_ids[!(ensembl_ids %in% valid_keys)]

# Print invalid keys, if any
print(invalid_keys)

# Now, your DESeq result should contain a new column with official gene symbols
gene_symbols <- select(org.Hs.eg.db, keys = ensembl_ids, columns = c("SYMBOL"), keytype = "ENSEMBL")

# Replace the Ensemble IDs in your DESeq results with gene symbols
deseq_result$gene_symbols_column_name <- gene_symbols 
