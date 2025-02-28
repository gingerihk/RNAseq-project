library(readr)
library(rtracklayer)
library(dplyr)
library(R.utils)
library(WGCNA)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(GEOquery)
library(tidyr)
library(genefilter)
library(ggrepel)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(pheatmap)
library(DescTools)
library(tinytex)

# Import table 
df2 <-readr::read_tsv("/Users/andreeaiuhaniak/Desktop/Software-Project/table_counts.tsv") 

# Download GTF file
gtf_url <- "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
dest_dir <- "/Users/andreeaiuhaniak/Desktop/Software-Project/"
dest_file <- paste0(dest_dir, "Homo_sapiens.GRCh38.110.gtf.gz")
download.file(gtf_url, destfile = dest_file, mode = "wb")
gunzip(dest_file, overwrite = TRUE)

# Load GTF file
gtf_data <- import("/Users/andreeaiuhaniak/Desktop/Software-Project/Homo_sapiens.GRCh38.110.gtf")
gtf_df <- as.data.frame(gtf_data)


# Extract column from GTF file 
transcript_to_gene <- gtf_data %>%
  as.data.frame() %>%
  filter(type == "transcript") %>%
  select(seqnames, transcript_id, gene_name, gene_id) %>%
  distinct()
colnames(transcript_to_gene) <- c("chromosome", "transcript_id", "gene_name", "gene_id")


# Merge counts table with transcript to genes 
  # adjust format of Transcript ID values
df2 <- df2 %>%
  mutate(Transcript_ID = sub("\\..*", "", Transcript_ID)) #target_id
df2$Transcript_ID <- as.factor(df2$Transcript_ID)
  # join df2 with transcript_to_gene
df2 <- left_join(df2, transcript_to_gene, by= c("Transcript_ID" = "transcript_id"))
  # sum transcript counts to get gene counts 
gene_counts_name <- df2 %>% group_by (gene_name) %>% summarise(across(where(is.numeric), sum, na.rm = TRUE)) 


# Meta_data table 

gse <- getGEO(filename = "/Users/andreeaiuhaniak/Downloads/GSE138746_series_matrix.txt")
meta_data_table <- data.frame(Sample_ID = gse$geo_accession, Sex = gse$`Sex:ch1`, Age = gse$`age:ch1`, DiseaseState = gse$`disease state:ch1`, Drug = gse$`drug:ch1`, Response = gse$`response:ch1`, CellType = gse$`cell type:ch1`)
meta_data_table <- readr::read_csv("/Users/andreeaiuhaniak/Desktop/meta_data_table.csv")

 # Create a vector of SRR numbers
srr_numbers <- paste0("SRR10260", 429:508) #sample number of monocytes 

 # Initialize a new column with NA
meta_data_table$SRR_ID <- NA

 # Assign SRR numbers only to rows where CellTypes == "monocytes"
meta_data_table$SRR_ID[meta_data_table$CellType == "Monocytes"] <- srr_numbers
meta_data_table_monocytes <- meta_data_table[81:160, ]

 #Adjust table 
meta_data_table_monocytes <- meta_data_table_monocytes[ ,-1] #removes first col Sample ID
meta_data_table_monocytes <- meta_data_table_monocytes[, c(ncol(meta_data_table_monocytes), 1:(ncol(meta_data_table_monocytes) - 1))] #moves last column SRR_ID to first position
meta_data_table_monocytes <- as.data.frame(meta_data_table_monocytes) #from tibble to df
row.names(meta_data_table_monocytes) <- meta_data_table_monocytes$SRR_ID #sets row names using SRR_IDs
meta_data_table_monocytes <- meta_data_table_monocytes[, -1] #removes SRR_ID col bc now it's stored as row names
meta_data_table_monocytes <- meta_data_table_monocytes[row.names(meta_data_table_monocytes) != "SRR10260461", ] #removes row with specific SRR ID

# 0 & 1 encoded meta_data table for correlation 

meta_table_monocytes_encoded <- meta_data_table_monocytes %>%
  mutate(
    Response = recode(Response, "no" = 1, "moderate" = 2, "good" = 3),
    Drug = recode(Drug, "Etanercept" = 1, "Adalimumab" = 2),
    Sex = recode(Sex, "f" = 1, "m" = 2)
  )

meta_table_monocytes_encoded <- meta_table_monocytes_encoded %>%
  select(-c(CellType, DiseaseState))


# matrix gene_name 

gene_counts_name <- as.data.frame(gene_counts_name) #so you can set gene_name as row name
gene_counts_name <- gene_counts_name[-nrow(gene_counts_name), ] # deletes last row 
row.names(gene_counts_name) <- gene_counts_name$gene_name # sets gene_name column as row names 
gene_counts_name <- gene_counts_name %>% select(-gene_name) #deletes column gene_name
gene_counts_name_m <- as.matrix(gene_counts_name) # deseq2 requires matrix 


#DESeq2 Gene Name 
dds_1 <- DESeqDataSetFromMatrix(round(gene_counts_name_m),   #nr rows in meta = columns df2_m 
                              meta_data_table_monocytes,
                              design = ~Drug + Response)
dds_1 <- DESeq(dds_1) #is variance stabilization necessary?
norm_counts_1 <- counts(dds_1, normalized = TRUE)
norm_counts_1_t <- t(norm_counts_1)

# chose soft threshold power 
powers <- c(1:20)  # Range of soft-thresholding powers
sft <- pickSoftThreshold(networkType= "signed", norm_counts_1_t, powerVector=powers, verbose=5)
#write.csv(gene_counts_name_t, file = "~/Desktop/gene_counts_name_wgcna.csv", row.names = FALSE)
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

picked_power <- 9
netwk <- blockwiseModules(norm_counts_1_t,  
                          power = picked_power,   
                          networkType = "signed",  
                          deepSplit = 3,         
                          corType = "pearson",
                          maxBlockSize = 5000,
                          minModuleSize = 30,
                          mergeCutHeight = 0.25,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "PW_9_30k_TOM",
                          numericLabels = TRUE,
                          verbose = 3)
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors) #only colors 

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

# Relate modules to samples 
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors) #gene id + color 
)

  # How many genes in each module 
module_counts <- table(module_df$colors) 
module_counts_df <- as.data.frame(module_counts)
colnames(module_counts_df) <- c("Module", "Gene_Count")


# Merging of modules whose expression profiles are very similar 

MEList <- moduleEigengenes(norm_counts_1_t, colors = netwk$colors, impute = TRUE)
MEs <- MEList$eigengenes
MEDissThres = 0.9  # modules with >= 0.1 correlation coefficient will be merged into one module
# Modules with a correlation ≥ (1 - MEDissThres) are merged 


merge = mergeCloseModules(norm_counts_1_t, netwk$colors, cutHeight = MEDissThres, verbose = 3)  # identifies modules that are highly correlated and combines them into a single module
# merge is a list containing updated module assignments and eigengenes after merging
mergedColors_2 <- labels2colors(merge$colors) #only colors 

# cluster dendrogram for new merged modules
plotDendroAndColors(netwk$dendrograms[[1]], cbind(mergedColors[netwk$blockGenes[[1]]], mergedColors_2[netwk$blockGenes[[1]]]),
                    c("Module Colors", "Merged Modules"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

  # relate new modules to sample + gene number in each module 
module_df_merged <- data.frame(
  gene_id = names(merge$colors),
  colors = labels2colors(merge$colors) # base for enrichment analysis 
)
module_counts_new <- table(module_df_merged$colors)
module_counts_new_df <- as.data.frame(module_counts_new)
colnames(module_counts_new_df) <- c("Module", "Gene_Count")
 
#save df as csv
write.csv(module_df_merged, file = "~/Desktop/module_df_reduced.csv", row.names = FALSE)


# bar plot all modules 
ggplot(module_counts_new_df, aes(x=reorder(Module, -Gene_Count), y=Gene_Count, fill=Module)) +
  geom_bar(stat="identity", color="black") +  # Create bars with black border
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none") +  # Rotate x-axis labels for readability
  labs(title="Number of Genes per Module", x="Module", y="Gene Count") +
  scale_fill_identity()
# bar plot modules with <2000 genes
module_counts_new_df_1 <- module_counts_new_df[-c(14, 23, 26, 33, 34), ] #delete rows because of high value
ggplot(module_counts_new_df_1, aes(x=reorder(Module, -Gene_Count), y=Gene_Count, fill=Module)) +
  geom_bar(stat="identity", color="black") +  # Create bars with black border
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none") +  # Rotate x-axis labels for readability
  labs(title="Number of Genes per Module", x="Module", y="Gene Count") +
  scale_fill_identity()


# MODULE-TRAIT 

nGenes <- ncol(norm_counts_1_t)
nSamples <- nrow(norm_counts_1_t)

# calculate eigengene values 
MEs0 <- moduleEigengenes(norm_counts_1_t, mergedColors_2)$eigengenes
MEs <- orderMEs(MEs0)


# Compute Pearson correlation between module eigengenes (MEs) and traits
modOXCor = cor(MEs, meta_table_monocytes_encoded[, colnames(meta_table_monocytes_encoded)[1:4]], use = "p")

# Format the correlation and p-values for display
textMatrix1 = paste(signif(modOXCor, 2))
dim(textMatrix1) = dim(modOXCor)

# Generate correlation heatmap (ComplexHeatmap)
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
Heatmap(modOXCor, name = "Correlation", 
        col = col_fun,
        column_names_side = "top",
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10),  # Adjust row label size
        column_names_gp = gpar(fontsize = 10),  # Adjust column label size
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", modOXCor[i, j]), x, y, gp = gpar(fontsize = 8))
        })

# Box plot eigengenes values grey based on sex 
MEs_grey <- MEs$MEgrey
eigengene_grey <- data.frame(
  Sample_ID <- rownames(meta_data_table_monocytes),
  Eigengene_value <- MEs_grey,
  Sex <- meta_data_table_monocytes$Sex
)

ggplot(eigengene_grey, aes(x = Sex, y = Eigengene_value, fill = Sex)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Grey Module Eigengene Values by Sex",
       x = "Sex", y = "Eigengene Value") +
  scale_fill_manual(values = c("m" = "blue", "f" = "pink"))

#Box = IQR (middle 50% data); Line = median; whiskers = min/max; dots = outliers 

# Box plot eigengenevalues red based on response
MEs_red <- MEs$MEred
eigengene_red <- data.frame(
  Sample_ID <- rownames(meta_data_table_monocytes),
  Eigengene_val <- MEs_red,
  Response <- meta_data_table_monocytes$Response
)

ggplot(eigengene_red, aes(x = Response, y = Eigengene_val, fill = Response)) +
  geom_boxplot() + 
  theme_minimal() +
  labs(title = "Red Module Eigengene Values by Response",
       x = "Response", y = "Eigengene Value") + 
  scale_fill_manual(values = c("no" = "blue", "good" = "magenta", "moderate" = "purple"))

#small eigenvalues -> low gene expression lvls; genes in the module do not vary much in expression across samples. 

# TranscriptID -> WGCNA
  # RED MODULE

red_genes <- module_df_merged %>% 
  filter(colors == "red") %>% 
  pull(gene_id)
print(red_genes)

red_genes_df <- data.frame(red_genes)
united_df <- red_genes_df %>%
  left_join(df2, by = c("red_genes" = "gene_name")) # merges df2 and red_genes_df based on genes present in red_genes_df
united_df <- united_df[, -((ncol(united_df) - 1):ncol(united_df))] # deletes last two columns(chromosome, gene_id)
row.names(united_df) <- united_df$Transcript_ID #sets row names by column 
united_df$Transcript_ID <- NULL #deletes Transcript ID
 united_df$red_genes <- NULL #deletes column red_genes (don't need it for wgcna)
united_matrix <- as.matrix(united_df) #?

   # normalization 
dds_2 <- DESeqDataSetFromMatrix(round(united_df),   #nr rows in meta = columns df2_m 
                                meta_data_table_monocytes,
                                design = ~Drug + Response)
dds_2 <- DESeq(dds_2) 
norm_counts_2 <- counts(dds_2, normalized = TRUE)
norm_counts_2_t <- t(norm_counts_2) #WGCNA input

  # soft threshold 
powers_1 <- c(1:20)  # Range of soft-thresholding powers
sft_1 <- pickSoftThreshold(networkType= "signed", norm_counts_2_t, powerVector=powers_1, verbose=5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft_1$fitIndices[, 1],
     -sign(sft_1$fitIndices[, 3]) * sft_1$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft_1$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft_1$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft_1$fitIndices[, 1],
     sft_1$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft_1$fitIndices[, 1],
     sft_1$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power <- 9
netwk_1 <- blockwiseModules(norm_counts_2_t,  
                          power = picked_power,   #deepsplit - higher. correlation, max blocksize smaller, 
                          networkType = "signed",   #normalized data for wgcna 
                          deepSplit = 4,          #main genes in each module; co splicing: use similar approach to do the co-splicing analysis
                          corType = "pearson",
                          maxBlockSize = 4000,
                          minModuleSize = 30,
                          mergeCutHeight = 0.25,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "PW_9_25k_TOM",
                          numericLabels = TRUE,
                          verbose = 3)
mergedcolors = labels2colors(netwk_1$colors)
plotDendroAndColors(
  netwk_1$dendrograms[[1]],
  mergedcolors[netwk_1$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

# Merging of modules whose expression profiles are very similar (red module)

MEList1 <- moduleEigengenes(norm_counts_1_t, colors = netwk$colors, impute = TRUE)
MEs1 <- MEList1$eigengenes
MEDiss1 = 1-cor(MEs1)
METree1 = hclust(as.dist(MEDiss1), method = "average");
MEDissThres1 = 0.9  # modules with >= 0.1 correlation coefficient will be merged into one module
# Modules with a correlation ≥ (1 - MEDissThres) are merged 
# lower diss thres -> fewer modules merged

merge1 = mergeCloseModules(norm_counts_2_t, netwk_1$colors, cutHeight = MEDissThres1, verbose = 3) #  identifies modules that are highly correlated and combines them into a single module
#merge is a list containing updated module assignments and eigengenes after merging
#netwk$colors has gene names!


mergedColors_red <- merge1$colors  # New merged module assignments
mergedColors_red1 <- labels2colors(merge1$colors) #only colors 
mergedColors_red_df <- as.data.frame(mergedColors_red1) # gene and number as color rn/ just colors ; also no need 
mergedMEs <- merge$newMEs # New eigengenes for merged modules (no need bc you dont use it)

plotDendroAndColors(netwk_1$dendrograms[[1]], cbind(mergedcolors[netwk_1$blockGenes[[1]]], mergedColors_red1[netwk_1$blockGenes[[1]]]),
                    c("Module Colors", "Merged Modules"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# Relate new modules to sample + gene number in each module 
red_module_df_merged <- data.frame(
  gene_id = names(merge1$colors),
  colors = labels2colors(merge1$colors) #base for enrichment analysis 
)
red_module_counts_new <- table(red_module_df_merged$colors)
red_module_counts_new_df <- as.data.frame(red_module_counts_new)
colnames(red_module_counts_new_df) <- c("Module", "Gene_Count")

# MODULE-TRAIT TRANSCRIPT ID RED MODULE

nGenes1 <- ncol(norm_counts_2_t)
nSamples1 <- nrow(norm_counts_2_t)

# calculate eigengene values 
MEs0_ <- moduleEigengenes(norm_counts_2_t, mergedColors_red1)$eigengenes
MEs_ <- orderMEs(MEs0_)

# Compute Pearson correlation between module eigengenes (MEs) and traits
modOXCor1 = cor(MEs_, meta_table_monocytes_encoded[, colnames(meta_table_monocytes_encoded)[1:4]], use = "p")

# Compute Pearson correlation p-values
modOXP1 = corPvalueStudent(modOXCor1, nSamples1)

# Format the correlation and p-values for display
#textMatrix2 = paste(signif(modOXCor1, 2), "\n(", signif(modOXP1, 1), ")", sep = "")
#dim(textMatrix2) = dim(modOXCor1)

# Convert matrices to data frames for easier handling
modOXCor_results1 <- as.data.frame(modOXCor1)
modOXP_results1 <- as.data.frame(modOXP1)

# Find the minimum correlation value
min(modOXCor_results1)

# Generate correlation heatmap for red module (ComplexHeatmap)
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
Heatmap(modOXCor1, name = "Correlation", 
        col = col_fun,
        column_names_side = "top",
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10),  # Adjust row label size
        column_names_gp = gpar(fontsize = 10),  # Adjust column label size
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", modOXCor1[i, j]), x, y, gp = gpar(fontsize = 8))
        })

# bar plot all submodules in red modules after WGCNA and merging 
ggplot(red_module_counts_new_df, aes(x=reorder(Module, -Gene_Count), y=Gene_Count, fill=Module)) +
  geom_bar(stat="identity", color="black") +  # Create bars with black border
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none", plot.title = element_text(hjust = 0.5)) +  # Rotate x-axis labels for readability
  labs(title="Number of Transcripts per Module", x="Module", y="Transcript Counts") +
  scale_fill_identity()

# bar plot modules with <2000 genes
red_module_counts_new_df_1 <- red_module_counts_new_df[-c(16, 30, 26, 12, 37), ] #delete rows because of high value
ggplot(red_module_counts_new_df_1, aes(x=reorder(Module, -Gene_Count), y=Gene_Count, fill=Module)) +
  geom_bar(stat="identity", color="black") +  # Create bars with black border
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none", plot.title = element_text(hjust = 0.5)) +  # Rotate x-axis labels for readability
  labs(title="Number of Transcripts per Module", x="Module", y="Transcript Counts") +
  scale_fill_identity()



# Splicing analysis
 # Splicing factors download and import 

url <- "http://srv00.recas.ba.infn.it/SpliceAidF/download.php?t=f&s=all"
destination <- "/Users/andreeaiuhaniak/Desktop/Software-Project/"
destination_file <- paste0(destination, "spliceaidf_data.tsv")
download.file(url, destfile = destination_file, mode = "w")

splice_data <- read_tsv("/Users/andreeaiuhaniak/Desktop/Software-Project/spliceaidf_data.tsv")
splice_data_df <- as.data.frame(splice_data)
# matching_genes <- subset(splice_data_df, Gene %in% genes_df$combined_module_genes)  # checks if there are matching splice factors in chosen modules
matching_genes_all_modules <- subset(splice_data_df, Gene %in% rownames(norm_counts_1_df))  # matching splice factors in all modules 

# matching_splicing_genes <- as.data.frame(splice_data_df[, 2]) #only splicing factors in df 
norm_counts_1_df <- as.data.frame(norm_counts_1)
norm_counts_3 <- rownames_to_column(norm_counts_1_df, var = "Gene_name") # row names to column "gene_name"
norm_counts_3 <- subset(norm_counts_3, Gene_name %in% matching_genes_all_modules$Gene)
which_modules <- subset(module_df_merged, gene_id %in% matching_genes_all_modules$Gene)

common_genes <- intersect(norm_counts_3$Gene_name, splice_data_df$Gene) #common genes between expression data and splicing factors
print(common_genes)

#WGCNA splicing factors 

row.names(norm_counts_3) <- norm_counts_3$Gene_name
norm_counts_3 <- norm_counts_3 %>% select(-Gene_name)
norm_counts_3_t <- t(norm_counts_3) # genes are columns 
picked_power <- 9
netwk_2 <- blockwiseModules(norm_counts_3_t,  
                            power = picked_power,   #deepsplit - higher. correlation, max blocksize smaller, 
                            networkType = "signed",   #normalized data for wgcna 
                            deepSplit = 3,          #main genes in each module; co splicing: use similar approach to do the co-splicing analysis
                            corType = "pearson",
                            maxBlockSize = 5000,
                            minModuleSize = 10,
                            mergeCutHeight = 0.25,
                            saveTOMs = TRUE,
                            saveTOMFileBase = "PW_8_366_TOM",
                            numericLabels = TRUE,
                            verbose = 3).  # not possible, too little genes, only grey modules 
mergedcolors1 <- labels2colors(netwk_2$colors)
plotDendroAndColors(
  netwk_2$dendrograms[[1]],
  mergedcolors1[netwk_2$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

# Splicing factors and counts correlation + plot; useful for identifying genes whose expression correlates with phenotypic traits
  # Compute correlation matrix
cor_matrix <- cor((norm_counts_3_t), meta_table_monocytes_encoded, method = "pearson", use = "pairwise.complete.obs")
pheatmap(cor_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Correlation between splicing factors expression and traits")

# Splicing factors + colors plot
  # Number of genes per Module
ggplot(which_modules, aes(x = colors, fill = colors)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Number of splicing factors per Module",
       x = "Module Color",
       y = "Number of Genes") +
  scale_fill_identity() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

  # Gene-Module Association 
ggplot(which_modules, aes(x = gene_id, y = colors, fill = colors)) +
  geom_tile() +
  theme_minimal() +
  labs(title = "Splicing factors-Module Association",
       x = "Gene",
       y = "Module Color") +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none", plot.title = element_text(hjust = 0.5))


## Likelihood of transcript from same gene to be in same module
# how many transcripts for each gene in red module 
red_gene_transcripts <- df2 %>%.  # turns Gene ID into Transcript ID by summing 
filter(gene_name %in% red_genes) %>%  # keeps only genes present in red_genes
  group_by(gene_name) %>%     # prepares data for summarization 
  summarise(TranscriptCount = n()) %>%  #Counts the number of transcripts for each gene
  arrange(desc(TranscriptCount))

red_genes_df_1 <- data.frame(red_genes)
united_df_1 <- red_genes_df_1 %>%
  left_join(df2, by = c("red_genes" = "gene_name")) %>% # merges df2 and red_genes_df based on genes present in red_genes_df
  left_join(red_module_df_merged, by = c("Transcript_ID" = "gene_id")) # merges red_module_df_merged that has transcripts and color swith the rest
united_df_1 <- united_df_1 %>%
  select(1:2, last_col(), 3:(ncol(.)-1)) # move last column (colors) to third position 
united_df_1 <- united_df_1[, -((ncol(united_df_1) - 1):ncol(united_df_1))] #deletes last two columns (chr, geneid)
united_df_1 <- united_df_1 %>% 
  left_join(red_gene_transcripts, by = c("red_genes" = "gene_name"))
united_df_1 <- united_df_1 %>%
  select(1:2, last_col(), 3:(ncol(.)-1))
united_df_1_filtered <- united_df_1 %>% group_by(red_genes) %>% filter(n() > 1) #filters out genes with 1 transcript only

# Idenfity transcript pair per gene 
# Count pairs of transcripts per gene and module
same_module_pairs <- united_df_1_filtered %>%
  group_by(red_genes, colors) %>%
  summarise(PairsInSameModule = choose(n(), 2)) %>%
  ungroup() %>%
  summarise(TotalSameModulePairs = sum(PairsInSameModule, na.rm = TRUE))

# Calculate total transcript pairs per gene
total_pairs <- united_df_1_filtered %>%
  group_by(red_genes) %>%
  summarise(TotalPairs = choose(n(), 2)) %>%
  ungroup() %>%
  summarise(TotalPairsSum = sum(TotalPairs, na.rm = TRUE))

# Calculate the expected pairs in the same module
num_modules <- length(unique(united_df_1_filtered$colors))
expected_pairs <- total_pairs$TotalPairsSum / num_modules


# Contingency Table
contingency_table <- matrix(c(same_module_pairs$TotalSameModulePairs,
                              total_pairs$TotalPairsSum - same_module_pairs$TotalSameModulePairs,
                              expected_pairs,
                              total_pairs$TotalPairsSum - expected_pairs),
                            nrow = 2,
                            byrow = TRUE)
colnames(contingency_table) <- c("Same Module", "Different Module")
rownames(contingency_table) <- c("Observed", "Expected")

chisq_test <- chisq.test(contingency_table)
print(chisq_test) # very significant; reject H0 which means strong likelihood that transcripts that come from the same gene are in the same module 

# HGS: gene with most transcripts 
gene_most_transcript <- which.max(united_df_1_filtered$TranscriptCount) #row 9831
subset_united_df_filtered <- united_df_1_filtered[united_df_1_filtered$red_genes == "HGS", c("red_genes", "colors")]
transcript_number <- as.data.frame(table(subset_united_df_filtered$colors))

ggplot(transcript_number, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 45)) + #names on X axis tilted 
  labs(
    title = paste("Transcript Module Distribution for HGS"),
    x = "Module Color",
    y = "Transcript Number"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + #title in the middle
  scale_fill_identity()

#expression of genes in magenta
# filter HGS genes 
subset_united_df_filtered_expression <- united_df_1_filtered[united_df_1_filtered$red_genes == "HGS", ]
# filter rows for magenta transcripts & drop unwanted columns 
magenta_module <- subset_united_df_filtered_expression %>% 
  ungroup() %>%
  filter(colors == "magenta") %>% 
  select(-'red_genes', -'TranscriptCount', -'colors')
magenta_module_long <- magenta_module %>% #pivot data to long format
  pivot_longer(
    cols = starts_with("SRR"),
    names_to = "sample",
    values_to = "expression"
  )
# calculate aaverage expression
average_expression_magenta <- magenta_module_long %>% group_by(Transcript_ID) %>% summarise(avg_expression = mean(expression, na.rm = TRUE))

# plot the average expression for each transcript id
ggplot(average_expression_magenta, aes(x = Transcript_ID, y = avg_expression)) +
  geom_bar(stat = "identity", fill = "magenta") +
  labs(title = "Average Expression of Transcripts in Module Magenta",
       x = "Transcript ID",
       y = "Average Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#expression for genes in pink 
# Filter for pink module and drop unwanted columns
pink_module <- subset_united_df_filtered_expression %>% 
  ungroup() %>%
  filter(colors == "pink") %>% 
  select(-red_genes, -TranscriptCount, -colors)

# Pivot the data to long format
pink_module_long <- pink_module %>% 
  pivot_longer(
    cols = starts_with("SRR"),
    names_to = "sample",
    values_to = "expression"
  )

# Calculate average expression for each transcript
average_expression_pink <- pink_module_long %>% 
  group_by(Transcript_ID) %>% 
  summarise(avg_expression = mean(expression, na.rm = TRUE))

# Plot the average expression for each transcript
ggplot(average_expression_pink, aes(x = Transcript_ID, y = avg_expression)) +
  geom_bar(stat = "identity", fill = "pink") +
  labs(title = "Average Expression of Transcripts in Module Pink",
       x = "Transcript ID",
       y = "Average Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# expression for genes in midnightblue
# Filter for midnight blue module and drop unwanted columns
midnightblue_module <- subset_united_df_filtered_expression %>% 
  ungroup() %>%
  filter(colors == "midnightblue") %>% 
  select(-red_genes, -TranscriptCount, -colors)

# Pivot the data to long format
midnightblue_module_long <- midnightblue_module %>% 
  pivot_longer(
    cols = starts_with("SRR"),
    names_to = "sample",
    values_to = "expression"
  )

# Calculate average expression for each transcript
average_expression_midnight_blue <- midnightblue_module_long %>% 
  group_by(Transcript_ID) %>% 
  summarise(avg_expression = mean(expression, na.rm = TRUE))

# Plot the average expression for each transcript
ggplot(average_expression_midnight_blue, aes(x = Transcript_ID, y = avg_expression)) +
  geom_bar(stat = "identity", fill = "midnightblue") +
  labs(title = "Average Expression of Transcripts in Module Midnightblue",
       x = "Transcript ID",
       y = "Average Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






















