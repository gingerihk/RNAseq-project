library(readr)
library(rtracklayer)
library(dplyr)
library(R.utils)
library(MEGENA)
library(WGCNA)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(GEOquery)
library(tidyr)
library(genefilter)
library(ggrepel)
library(matrixStats)

df<-readr::read_tsv("/Users/andreeaiuhaniak/Desktop/merged_counts.tsv") #old counts with all files
df2<-readr::read_tsv("/Users/andreeaiuhaniak/Desktop/table_counts.tsv") #new counts
df2<- df2[-1, ]
df3 <-readr::read_tsv("/Users/andreeaiuhaniak/Desktop/table_tpm.tsv")
df3<-df3[-1, ]

#Download GTF file
gtf_url <- "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
dest_dir <- "/Users/andreeaiuhaniak/Desktop/Software-Project/"
dest_file <- paste0(dest_dir, "Homo_sapiens.GRCh38.110.gtf.gz")
download.file(gtf_url, destfile = dest_file, mode = "wb")
gunzip(dest_file, overwrite = TRUE)

# Load GTF file
gtf_data <- import("/Users/andreeaiuhaniak/Desktop/Software-Project/Homo_sapiens.GRCh38.110.gtf")
gtf_df <- as.data.frame(gtf_data)

#Load MANE file 
mane_file <- "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz"
dest_dir <- "/Users/andreeaiuhaniak/Desktop/Software-Project/"
dest_file <- paste0(dest_dir, "MANE.GRCh38.v1.4.summary.txt.gz")
download.file(mane_file, destfile = dest_file, mode = "wb")

# Extract transcript-to-gene mapping with gene names
transcript_to_gene <- gtf_data %>%
  as.data.frame() %>%
  filter(type == "transcript") %>%
  select(seqnames, transcript_id, gene_name, gene_id) %>%
  distinct()
colnames(transcript_to_gene) <- c("chromosome", "transcript_id", "gene_name", "gene_id")

# Merge counts table with transcript to genes
df2 <- df2 %>%
  mutate(Transcript_ID = sub("\\..*", "", Transcript_ID)) #target_id
#cpm_df2 <- t(t(df2[-1]) / colSums(df2[-1])) * 1e6
#normalized_df2 <- cbind(Transcript_ID = df2$Transcript_ID, cpm_df2)

#gene id
df2 <- left_join(df2, transcript_to_gene, by= c("Transcript_ID" = "transcript_id"))
gene_counts <- df2 %>% group_by (gene_id) %>% summarise(across(where(is.numeric), sum, na.rm = TRUE))
gene_counts_name <- df2 %>% group_by (gene_name) %>% summarise(across(where(is.numeric), sum, na.rm = TRUE))
#gene_counts_filtered <- gene_counts[rowSums(norm_counts > 10) > 0, ]
merged_counts <- df2 %>%
  left_join(transcript_to_gene, by = c("Transcript_ID" = "transcript_id"))
merged_counts <- merged_counts %>%
  select(Transcript_ID, gene_name, chromosome, gene_id, everything())

#isoforms
isoform_counts <- as.data.frame(table(merged_counts$gene_name))
colnames(isoform_counts) <- c("gene_name", "Num_Isoforms")
merged_data <- merge(merged_counts, isoform_counts, by="gene_name")
merged_data <- merged_data[, c(1:3, ncol(merged_data), 4:(ncol(merged_data)-1))]
ggplot(merged_data, aes(x=Num_Isoforms)) +
  geom_histogram(binwidth=1, fill="steelblue", color="black") +
  labs(title="Distribution of Isoforms per Gene",
       x="Number of Isoforms per Gene", y="Frequency") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
count_matrix_long <- melt(as.data.frame(df2), 
                          variable.name = "sample", 
                          value.name = "count", 
                          id.vars = "Transcript_ID") #df with transcript, sample, count
class(count_matrix_long)

#df3 <- df3 %>%
  #mutate(Transcript_ID = sub("\\..*", "", Transcript_ID)) #target_id

#merged_tpm <- df3 %>%
  #left_join(transcript_to_gene, by = c("Transcript_ID" = "transcript_id"))
#merged_tpm <- merged_tpm %>%
  #select(Transcript_ID, gene_name, chromosome, everything())


##split CD14

filtered_CD14 <- merged_counts %>%
  select(1:3, 46:106)

# Meta_data table 

gse <- getGEO(filename = "/Users/andreeaiuhaniak/Downloads/GSE138746_series_matrix.txt")
meta_data_table <- data.frame(Sample_ID = gse$geo_accession, Sex = gse$`Sex:ch1`, Age = gse$`age:ch1`, DiseaseState = gse$`disease state:ch1`, Drug = gse$`drug:ch1`, Response = gse$`response:ch1`, CellType = gse$`cell type:ch1`)

 # Create a vector of SRR numbers
srr_numbers <- paste0("SRR10260", 429:508)

 # Initialize a new column with NA
meta_data_table$SRR_ID <- NA

 # Assign SRR numbers only to rows where CellTypes == "monocytes"
meta_data_table$SRR_ID[meta_data_table$CellType == "Monocytes"] <- srr_numbers
meta_data_table_monocytes <- meta_data_table[81:160, ]

##Long format df2 (SRR rows, Transcript columns)

df2_long <- df2 %>%
  pivot_longer(
    cols = -Transcript_ID,  # All columns except `transcript_ID`
    names_to = "SRR_ID",    # Create a column `SRR_ID` from column names
    values_to = "counts"    # Create a column `counts` for values
  )

df2_merged <- df2_long %>%
  inner_join(meta_data_table_monocytes, by = "SRR_ID")

#Visualize how transcript counts vary by response categories
boxplot <- ggplot(df2_merged, aes(x = Response, y = counts)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.3, color = "blue") +
  theme_minimal() +
  labs(x = "Response", y = "Counts", title = "Transcript Counts by Response")



#DESeq2 normalization 

#change to matrix gene_ID
gene_counts <- as.data.frame(gene_counts) #so you can set transcript_id as rowname
gene_counts <- gene_counts[-nrow(gene_counts), ]
row.names(gene_counts) <- gene_counts$gene_id

gene_counts <- gene_counts %>% select(-gene_id) #deletes column Transcript_id
gene_counts_m <- as.matrix(gene_counts) 
gene_counts_m <- gene_counts_m[-nrow(gene_counts_m), ] #removes last row bc it has NA
#gene_counts_m_filtered <- gene_counts_m[rowSums(gene_counts_m > 0) > 0, ]

gene_counts_name <- as.data.frame(gene_counts_name) #so you can set transcript_id as rowname
gene_counts_name <- gene_counts_name[-nrow(gene_counts_name), ]
row.names(gene_counts_name) <- gene_counts_name$gene_name
gene_counts_name <- gene_counts_name %>% select(-gene_name) #deletes column Transcript_id
gene_counts_name_m <- as.matrix(gene_counts_name) 
gene_counts_name_m <- gene_counts_name_m[-nrow(gene_counts_name_m), ] #removes last row bc it has NA
#gene_counts_name_m_filtered <- gene_counts_name_m[rowSums(gene_counts_name_m > 0) > 0, ]
gene_counts_name_t = t(gene_counts_name) #input wgcna
#change to matrix

df2$Transcript_ID <- gsub("[.].*", "", df2$Transcript_ID) 
df2 <- as.data.frame(df2) #so you can set transcript_id as rowname
row.names(df2) <- df2$Transcript_ID  #sets names of rows 
df2 <- df2 %>% select(-Transcript_ID) #deletes column Transcript_id
df2_m <- as.matrix(df2) 
df2_m <- df2_m[-nrow(df2_m), ] #removes last row bc it has NA

meta_data_table_monocytes <- meta_data_table_monocytes[ ,-1]
meta_data_table_monocytes <- meta_data_table_monocytes[, c(ncol(meta_data_table_monocytes), 1:(ncol(meta_data_table_monocytes) - 1))]
row.names(meta_data_table_monocytes) <- meta_data_table_monocytes$SRR_ID
meta_data_table_monocytes <- meta_data_table_monocytes[, -1]
meta_data_table_monocytes <- meta_data_table_monocytes[row.names(meta_data_table_monocytes) != "SRR10260461", ]

#DESeq2 Gene ID
dds_1 <- DESeqDataSetFromMatrix(round(gene_counts_name_m_filtered),   #nr rows in meta = columns df2_m 
                              meta_data_table_monocytes,
                              design = ~Drug + Response)
dds_1 <- DESeq(dds_1)
norm_counts_1 <- counts(dds_1, normalized = TRUE)
#scaled_norm <- t(scale(norm_counts_1))
results_1 <- resultsNames(dds_1)
vsd_1<- varianceStabilizingTransformation(dds_1) #stabilizes the variance across the range of counts; imp for extreme expression genes; useful for PCA and clustering
wpn_vsd_1 <- getVarianceStabilizedData(dds_1)
rv_wpn_1 <- rowVars(wpn_vsd_1) #calculates variance for each gene
#filtered_genes_1 <- vsd_1[rv_wpn_1 > 2, ] #filter genes above 2
filtered_genes_1
#for WGCNA
stab_variance <- varianceStabilizingTransformation(dds_1)
wgcna_input <- assay(stab_variance) #normalized expression values
gene_variance <- rowVars(wgcna_input)
var_threshold <- quantile(gene_variance, 0.70)
filtered_wgcna_input <- wgcna_input[gene_variance > var_threshold, ] #keeps only informative genes
filtered_wgcna_input_t <- t(filtered_wgcna_input)
powers <- c(1:20)  # Range of soft-thresholding powers
sft <- pickSoftThreshold(networkType= "signed", gene_counts_name_t, powerVector=powers, verbose=5)
write.csv(gene_counts_name_t, file = "~/Desktop/gene_counts_name_wgcna.csv", row.names = FALSE)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")
# Plot the scale-free topology fit index
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale-Free Topology Fit", type="b",
     main="Soft Threshold Selection")
picked_power <- 9
netwk <- blockwiseModules(gene_counts_name_t,  
                          power = picked_power, 
                          networkType = "signed",
                          deepSplit = 2,
                          minModuleSize = 30,
                          mergeCutHeight = 0.25,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "PW_9_30k_TOM",
                          numericLabels = TRUE,
                          verbose = 3)
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

#Relate modules to samples 
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors) #gene id + color 
)

MEs0 <- moduleEigengenes(gene_counts_name_t, mergedColors)$eigengenes
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
MEs0$Sample <- row.names(MEs0) #adds new column called Sample with row names

mME = MEs0 %>%
  pivot_longer(-Sample) %>%
  mutate(
    name = gsub("ME", "", name),  # Remove "ME" from module names
    name = factor(name, levels = module_order))
mME %>% ggplot(., aes(x=Sample, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "green",
    high = "purple",
    mid = "gray",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
#Examine expression profiles 
  #pick modules of interest
modules_of_interest = c("tan", "red", "blue")

  #pull out genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)
row.names(module_df) = module_df$gene_id
subexpr <- filtered_wgcna_input_t[submod$gene_id, ] #contains expression data for selected module genes

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )
module_colors <- c(
  "blue" = "blue",
  "red" = "red",
  "tan" = "tan"
)

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  scale_color_manual(values = module_colors) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "samples",
       y = "normalized expression")

#PCA
pca_results <- prcomp(scaled_norm, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_results$x[, 1:2])
pca_df$Response <- meta_data_table_monocytes$Response
pca_df$Drug <- meta_data_table_monocytes$Drug
pca_df_A <- pca_df %>% filter(Drug == "Adalimumab")
pca_df_E <- pca_df %>% filter(Drug == "Etanercept")
colnames(pca_df) <- c("PC1", "PC2")

head(pca_df)

ggplot(pca_df_A, aes(x = PC1, y = PC2, color = Response, shape = Drug)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "PCA Analysis", x = "PC1", y = "PC2")
ggplot(pca_df_E, aes(x = PC1, y = PC2, color = Response, shape = Drug)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "PCA Analysis", x = "PC1", y = "PC2")

res_good_vs_no <- results(dds_1, contrast = c("Response", "good", "no"))
res_good_vs_no <- res_good_vs_no %>%
  left_join(meta_data_table_monocytes, by = "gene_name") #or gene_id
volcano_data <- as.data.frame(res_good_vs_no)
summary(res_good_vs_no)
# Apply thresholds for significance and fold change
volcano_data$significant <- ifelse(
  abs(volcano_data$log2FoldChange) > 1.5 & volcano_data$padj < 0.05, "Yes", "No"
)

# Color upregulated genes (log2FoldChange > 1.5) in blue, downregulated (log2FoldChange < -1.5) in red, 
# and non-significant in grey
volcano_data$color <- ifelse(
  volcano_data$significant == "Yes" & volcano_data$log2FoldChange > 1.5, "red",   # Red for upregulated genes
  ifelse(volcano_data$significant == "Yes" & volcano_data$log2FoldChange < -1.5, "blue", "grey")  # blue for downregulated, grey for non-significant
)
top_upregulated <- volcano_data[volcano_data$significant == "Yes", ] %>% 
  dplyr::arrange(desc(log2FoldChange)) %>% 
  head(5)

top_downregulated <- volcano_data[volcano_data$significant == "Yes", ] %>% 
  dplyr::arrange(log2FoldChange) %>% 
  head(5)
top_genes <- rbind(top_upregulated, top_downregulated)  # Combine top genes
 # Select the top 10 upregulated and top 10 downregulated genes
  top_upregulated <- volcano_data[volcano_data$significant == "Yes", ] %>% 
  dplyr::arrange(desc(log2FoldChange)) %>% 
  head(5)

top_downregulated <- volcano_data[volcano_data$significant == "Yes", ] %>% 
  dplyr::arrange(log2FoldChange) %>% 
  head(5)

top_genes <- rbind(top_upregulated, top_downregulated)  # Combine top 10 up/downregulated genes
# Create the volcano plot with annotations
  ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
  geom_point(alpha = 0.6, size = 2) + # Plot points
  scale_color_identity() +  # Use color column directly
  theme_minimal() +  # Clean theme
  labs(
    title = "Volcano Plot: Good vs No Response",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value"
  ) +
  # Add threshold lines
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "black", size = 1) +  # Fold change thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1) +  # P-value threshold
  # Annotate top 10 up and downregulated genes
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), 
                  size = 4, box.padding = 0.3, point.padding = 0.3, 
                  max.overlaps = 10) +  # Prevent excessive overlaps
  theme(
    plot.title = element_text(hjust = 0.5), 
    axis.title = element_text(size = 12), 
    axis.text = element_text(size = 10)
  )
# Create volcano plot
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
  geom_point(alpha = 0.6, size = 2) + # Plot points
  scale_color_identity(guide = "legend", labels = c(
    "blue" = "Downregulated genes",
    "red" = "Upregulated genes",
    "grey" = "Non-significant genes"
  )) + 
  theme_minimal() +  # Clean theme
  labs(
    title = "Volcano Plot: Good vs No Response",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), 
    axis.title = element_text(size = 12), 
    axis.text = element_text(size = 10)
  ) + geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "black", size = 1) +  # Fold change thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1) +  # P-value threshold
  # Annotate top 10 up and downregulated genes
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), 
                  size = 4, box.padding = 0.3, point.padding = 0.3, 
                  max.overlaps = 5) +  # Prevent excessive overlaps
  theme(
    plot.title = element_text(hjust = 0.5), 
    axis.title = element_text(size = 12), 
    axis.text = element_text(size = 10)
  )




# Set thresholds for significance and fold change
log2fc_threshold <- 1.5  # Adjust this threshold as needed
padj_threshold <- 0.05

# Apply thresholds to DESeq2 results
res_good_vs_no$significant <- ifelse(
  abs(res_good_vs_no$log2FoldChange) > log2fc_threshold & 
    res_good_vs_no$padj < padj_threshold, 
  "Yes", "No"
)

# View the updated results with significance column
head(res_good_vs_no)


# Filter significant genes
significant_genes_good_vs_no <- subset(
  res_good_vs_no, 
  padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold
)

# View significant genes
head(significant_genes_good_vs_no)

# Sort the significant genes by log2FoldChange in descending order
sorted_results_good_vs_no <- significant_genes_good_vs_no[order(-significant_genes_good_vs_no$log2FoldChange),]

# View sorted results
head(sorted_results_good_vs_no)

#ModeratevsNo
res_moderate_vs_no <- results(dds_1, contrast = c("Response", "moderate", "no"))
res_moderate_vs_no$significant <- ifelse(
  abs(res_moderate_vs_no$log2FoldChange) > log2fc_threshold & 
    res_moderate_vs_no$padj < padj_threshold, 
  "Yes", "No"
)

# View the updated results with significance column
head(res_moderate_vs_no)

# Filter significant genes
significant_genes_moderate_vs_no <- subset(
  res_moderate_vs_no, 
  padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold
)

# View significant genes
head(significant_genes_moderate_vs_no)

# Sort the significant genes by log2FoldChange in descending order
sorted_results_moderate_vs_no <- significant_genes_moderate_vs_no[order(-significant_genes_moderate_vs_no$log2FoldChange),]

# View sorted results
head(sorted_results_moderate_vs_no)




#DESeq2 Transcript ID
dds <- DESeqDataSetFromMatrix(round(df2_m),   #nr rows in meta = columns df2_m 
                              meta_data_table_monocytes,
                              design = ~Drug + Response)
dds <- DESeq(dds) #normalized data
res <- resultsNames(dds)
norm_counts <- counts(dds, normalized = TRUE)
#res <- results(dds, name="Response_no_vs_good")  
results <- results(dds)
vsd <- varianceStabilizingTransformation(dds) #stabilizes the variance across the range of counts; imp for extreme expression genes; useful for PCA and clustering
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd) #calculates variance for each gene
filtered_genes <- vsd[ rv_wpn > 2, ] #filter genes above 2
plotPCA(filtered_genes, intgroup=c("Drug", "Response"))
hist(rv_wpn, breaks=50, main="Variance Distribution", col="skyblue")
#summary(rv_wpn)
q80_wpn <- quantile( rowVars(filtered_genes), .80) #80th percentile, keeps only the genes in the top 20% most variable 
expr_normalized <- wpn_vsd[ rv_wpn > q80_wpn, ]
plotPCA(filtered_genes, intgroup=c("Drug", "Response"))



