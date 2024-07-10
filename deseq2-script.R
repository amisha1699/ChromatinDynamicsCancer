# Load DESeq2 library
library(DESeq2)

# Read count data for all samples
wt_samples <- c("SRR8435995", "SRR8435996", "SRR8435997",
                "SRR19159298", "SRR19159299", "SRR19159300",
                "SRR24572364", "SRR24572365")

mut_samples <- c("SRR8435992", "SRR8435993", "SRR8435994",
                 "SRR22192978", "SRR22192979", "SRR22192981",
                 "SRR22729521", "SRR22729522", "SRR22729523",
                 "SRR24442519", "SRR24442520", "SRR24442521")

# Read count data for all samples
wt_counts <- lapply(wt_samples, function(sample_id) {
  read.delim(paste0(sample_id, ".csv"), row.names = 1)
})

mut_counts <- lapply(mut_samples, function(sample_id) {
  read.delim(paste0(sample_id, ".csv"), row.names = 1)
})

# Combine replicates for each sample
wt_combined <- Reduce("+", wt_counts)
mut_combined <- Reduce("+", mut_counts)

# Create sample metadata
wt_replicates <- c(3, 3, 2)  # Number of replicates for each wild type sample
mut_replicates <- rep(3, 4)   # Number of replicates for each mutant sample

# Create sample metadata
sample_metadata <- data.frame(
  sampleName = c(rep(wt_samples[1], 3), rep(wt_samples[2], 3), rep(wt_samples[3], 2), 
                 rep(mut_samples[1], 3), rep(mut_samples[2], 3), rep(mut_samples[3], 3), rep(mut_samples[4], 3)),
  condition = c(rep("WT", 8), rep("Mutant", 12)),
  replicate = rep(rep(1:3, each = 3), times = c(2, 2, 1, 3, 3, 3, 3))  # Adjust the number of replicates
)



# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cbind(wt_combined, mut_combined),
                              colData = sample_metadata,
                              design = ~ condition + replicate)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract differential expression results
res <- results(dds)

# Perform downstream analysis and visualization as needed
