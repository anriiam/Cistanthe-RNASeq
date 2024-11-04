rm(list=ls())

#Commented out because I've already installed the program
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#  BiocManager::install(version = "3.19")
#BiocManager::install("DESeq2")
#BiocManager::install("maSigPro")
#BiocManager::install("DEGreport")


################ RNA DGE analysis Pipeline ################
# Anri Chomentowska June 2024

### Load required libraries
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(ashr)
library(ggpubr)
library(ggrepel)
library(stringr)
library(org.At.tair.db)
library(clusterProfiler)
library(GO.db)
library(AnnotationDbi)
library(DOSE)
library(enrichplot)
#library(biomaRt)
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)


######## Data importing and cleaning ########
#### Set working directory
# Where your outputs from FeatureCounts are downloaded on your local machine or Cluster
setwd("~/Dropbox/Yale/Research/RNAseq")


#### Load your data
mycounts <- read.table("finalcounts_data_stranded_geneid.txt", header = TRUE, as.is=T, skip = 1, row.names = 1)
mycounts_raw <- read.table("finalcounts_data_stranded_geneid.txt", header = TRUE, as.is=T, skip = 1, row.names = 1)
metadata <- read.csv("final_metadata.csv")

class(mycounts) # make sure these are data.frame
class(mycounts_raw)
class(metadata)

View(mycounts) # take a look at your data
View(mycounts_raw)
View(metadata)

# Clean up the count data
mycounts <- mycounts[ ,c(-1:-5)]
colnames(mycounts) <- gsub("Aligned.sortedByCoord.out.bam", "", colnames(mycounts), fixed = T)
colnames(mycounts) <- gsub("..star_results_new.", "", colnames(mycounts), fixed = T)
View(mycounts)

mycounts_raw <- mycounts_raw[ ,c(-1:-4)]
colnames(mycounts_raw) <- gsub("Aligned.sortedByCoord.out.bam", "", colnames(mycounts_raw), fixed = T)
colnames(mycounts_raw) <- gsub("..star_results_new.", "", colnames(mycounts_raw), fixed = T)
View(mycounts_raw)


## Subset Experiment 2: Pre_Flowering pre- and post-reproduction experiment
meta_wr <- metadata %>%
  dplyr::filter(Daylength == "Long")

myco_wr <- mycounts %>%
  dplyr::select(L9B6PM,
                L9B2PM,
                L9B10PM,
                L9B10AM,
                L8K6AM,
                L8K2PM,
                L8K2AM,
                L8K10AM,
                L8A6PM,
                L8A6AM,
                L8A2PM,
                L8A10PM,
                L8A10AM,
                L5I6PM,
                L5I6AM,
                L5I2PM,
                L5I2AM,
                L5I10PM,
                L5I10AM,
                L4L6PM,
                L4L6AM,
                L4L2PM,
                L4L2AM,
                L4L10PM,
                L4L10AM,
                L3I6PM,
                L3I6AM,
                L3I2PM,
                L3I2AM,
                L3I10PM,
                L3I10AM,
                L2Z6PM,
                L2Z6AM,
                L2Z2PM,
                L2Z2AM,
                L2Z10PM,
                L2Z10AM,
                L10O6PM,
                L10O6AM,
                L10O2PM,
                L10O2AM,
                L10O10PM,
                L10O10AM)

# Add a "Category" column
meta_wr <- meta_wr %>%
  mutate(Category = case_when(
    State %in% c("Vegetative", "Reproductive") ~ "Pre_Flowering",
    State %in% c("Flowering", "Fruiting") ~ "Post_Flowering",
    TRUE ~ as.character(State)
  ))

View(meta_wr)
View(myco_wr)

# Make sure the IDs are identical
names(myco_wr)
meta_wr$SampleID

names(myco_wr)==meta_wr$SampleID
all(names(myco_wr)==meta_wr$SampleID)

# Set correct order for time course and how data will be shown
meta_wr$Time <- factor(meta_wr$Time, levels = c("10AM","2PM","6PM","10PM","2AM","6AM"))
meta_wr$Category <- factor(meta_wr$Category, levels = c("Pre_Flowering","Post_Flowering"))

View(meta_wr)


#### Experiment 2: Pre_Flowering Pre- and Post-Reproductive Category ~ Time
#### Create DESeq2Dataset object
# Ignore the warning that comes up when it says "converting counts to integer mode..."
dds_wr <- DESeqDataSetFromMatrix(countData = round(myco_wr), 
                                 colData = meta_wr, 
                                 design = ~ Time + Category, 
                                 tidy = FALSE)
dds_wr # Take a look at the DESeq2 data structure

#### QC steps
# Perform the median of ratios method of normalization
# (This is also automatically performed by the DESeq() function, but you can take this output for use with other packages)
dds_wr <- estimateSizeFactors(dds_wr) # calculate normalization factor per sample (if significantly off from 1, then this means)

sizeFactors(dds_wr) #if very off from 1, then this means extreme outliers may exist. It doesn't here.

# retrieve the normalized counts matrix
normalized_counts_wr <- counts(dds_wr, normalized=TRUE)
write.table(normalized_counts_wr, file="normalized_counts_wr_NEW.txt", sep="\t", quote=F, col.names=NA) # Use this for analyses with other packages!
write.csv(meta_wr, file="metadata_wr_NEW.csv")

## Regularized log transformation (rlog) and Variance stabilizing transformation (vst) for data visualization
rld_wr <- rlog(dds_wr, blind=TRUE)
vst_wr <- vst(dds_wr, blind=TRUE)

# Extract the rlog and vst matrix from the object
rld_wr_mat <- assay(rld_wr)
vst_wr_mat <- assay(vst_wr)

# Compute pairwise correlation values
rld_wr_cor <- cor(rld_wr_mat)    ## cor() is a base R function
head(rld_wr_cor)   ## check the output of cor(), make note of the rownames and colnames

vst_wr_cor <- cor(vst_wr_mat)    ## cor() is a base R function
head(vst_wr_cor)   ## check the output of cor(), make note of the rownames and colnames

# Make first colum of metadata row names
meta_wr2 <- meta_wr[,-1]
rownames(meta_wr2) <- meta_wr[,1]
head(meta_wr2)

## Plot heatmap
pheatmap(rld_wr_cor, col = hmcol)
pheatmap(vst_wr_cor)

# Or, select the top 20 genes with the highest mean expression for heatmap
select_wr <- order(rowMeans(counts(dds_wr,normalized=TRUE)), decreasing=TRUE)[1:20]
anno_wr <- as.data.frame(colData(dds_wr)[,c("Time","Category")])

anno_wr$Time <- factor(anno_wr$Time, levels = c("10AM","2PM","6PM","10PM","2AM","6AM"))
anno_wr$Category <- factor(anno_wr$Category, levels = c("Pre_Flowering", "Post_Flowering"))

# Order the columns in the matrix based on Category and then Time
ordered_cols <- order(anno_wr$Category, anno_wr$Time)
anno_wr_ordered <- anno_wr[ordered_cols, ]

rld_wr_mat_ordered <- rld_wr_mat[, ordered_cols]
vst_wr_mat_ordered <- vst_wr_mat[, ordered_cols]

# Plot the heatmap of for the selected top 20 genes to visualize
pheatmap(rld_wr_mat_ordered[select_wr,], cluster_rows=FALSE, col = hmcol,
         cluster_cols=FALSE, annotation_col=anno_wr_ordered)

pheatmap(vst_wr_mat_ordered[select_wr,], cluster_rows=FALSE,
         cluster_cols=FALSE, annotation_col=anno_wr_ordered)

## Now, Plot a PCA 
plotPCA(rld_wr, intgroup="Category")
plotPCA(vst_wr, intgroup="Category")

# Different ways to explore more about the PCs (https://tavareshugo.github.io/data-carpentry-rnaseq/03_rnaseq_pca.html)
pca_wr <- prcomp(t(rld_wr_mat))
#pca_wr <- prcomp(t(vst_wr_mat)) #if you want to use vst
View(pca_wr)
df_wr <- cbind(meta_wr, pca_wr$x) # You can explore many ways to do this

pc_loadings_wr <- pca_wr$rotation
pc_loadings_wr <- pc_loadings_wr %>% 
  as_tibble(rownames = "gene")
View(pc_loadings_wr)

top_genes_wr <- pc_loadings_wr %>% 
  # select only the PCs we are interested in
  dplyr::select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows from PC1 and PC2 each
  top_n(10, loading) %>%
  # pull the gene column as a vector
  pull(gene)

top_genes_wr #list of the top genes (based on PC loadings)

# Filter the top loading genes' PC matrix and visualize the direction and magnitude of the eigenvalues
top_loadings_wr <- pc_loadings_wr %>% 
  filter(gene %in% top_genes_wr)

ggplot(top_loadings_wr) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "blue") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))

#another way to visualize and annotate PCA
pca_data_wr <- plotPCA(rld_wr, intgroup=c("Category","Time"), returnData=TRUE)
View(pca_data_wr)

percentVar_wr <- 100*round(attr(pca_data_wr, "percentVar"),2)               
makeLab <- function(x,pc) paste0("PC",pc,": ",x,"% variance")

ggplot(pca_data_wr) + geom_point(aes(x=PC1, y=PC2, color = Category, shape = Time), size=3) +
  xlab(makeLab(percentVar_wr[1],1)) + ylab(makeLab(percentVar_wr[2],2)) + 
  scale_color_manual(labels = c("Pre_Flowering", "Post_Flowering"), values=c("seagreen3", "magenta2")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_blank()) +
  ggtitle("PCA (rlog)")

#vst
pca_data_wr <- plotPCA(vst_wr, intgroup=c("Category","Time"), returnData=TRUE)
View(pca_data_wr)

percentVar_wr <- 100*round(attr(pca_data_wr, "percentVar"),2)               
makeLab <- function(x,pc) paste0("PC",pc,": ",x,"% variance")

ggplot(pca_data_wr) + geom_point(aes(x=PC1, y=PC2, color = Category, shape = Time), size=3) +
  xlab(makeLab(percentVar_wr[1],1)) + ylab(makeLab(percentVar_wr[2],2)) + 
  scale_color_manual(labels = c("Pre_Flowering", "Post_Flowering"), values=c("seagreen3", "magenta2")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_blank()) +
  ggtitle("PCA (vst)")


######## Differential gene expression (DEG) steps ########

## TEST 1 (testing just Category effects via LRT, for any differences over multiple time points)
# I think this is what I want to use!
dds_wr <- DESeqDataSetFromMatrix(countData = round(myco_wr), 
                                 colData = meta_wr, 
                                 design = ~ Time + Category, 
                                 tidy=FALSE)

dds_wr <- DESeq(dds_wr, test = "LRT", reduced = ~ Time)


# TEST 2:Treatment effects only 
#  ~ Time + Category + Time:Category
dds_wr <- DESeqDataSetFromMatrix(countData = round(myco_wr),
                                 colData = meta_wr, 
                                 design = ~ Category, 
                                 tidy=FALSE)

dds_wr <- DESeq(dds_wr, test = "LRT", reduced = ~1)
dds_wr <- DESeq(dds_wr) #Wald Test

dds_wr


# TEST 3.1 (p-value calculated from effect of Time, but between each Category): Post_Flowering
meta_po <- meta_wr %>%
  dplyr::filter(Category == "Post_Flowering")

myco_po <- myco_wr %>%
  dplyr::select(L9B6PM,
                L9B2PM,
                L9B10PM,
                L9B10AM,
                L5I6PM,
                L5I6AM,
                L5I2PM,
                L5I2AM,
                L5I10PM,
                L5I10AM,
                L4L6PM,
                L4L6AM,
                L4L2PM,
                L4L2AM,
                L4L10PM,
                L4L10AM,
                L2Z6PM,
                L2Z6AM,
                L2Z2PM,
                L2Z2AM,
                L2Z10PM,
                L2Z10AM)

dds_wr <- DESeqDataSetFromMatrix(countData = round(myco_po), 
                                 colData = meta_po, 
                                 design = ~ Time, 
                                 tidy=FALSE)
dds_wr <- DESeq(dds_wr, test = "LRT", reduced = ~1)

# TEST 3.2 (p-value calculated from effect of Time, but between each Category): Pre_Flowering
meta_pr <- meta_wr %>%
  dplyr::filter(Category == "Pre_Flowering")

myco_pr <- myco_wr %>%
  dplyr::select(L8K6AM,
                L8K2PM,
                L8K2AM,
                L8K10AM,
                L8A6PM,
                L8A6AM,
                L8A2PM,
                L8A10PM,
                L8A10AM,
                L3I6PM,
                L3I6AM,
                L3I2PM,
                L3I2AM,
                L3I10PM,
                L3I10AM,
                L10O6PM,
                L10O6AM,
                L10O2PM,
                L10O2AM,
                L10O10PM,
                L10O10AM)

dds_wr <- DESeqDataSetFromMatrix(countData = round(myco_pr), 
                                 colData = meta_pr, 
                                 design = ~ Time, 
                                 tidy=FALSE)
dds_wr <- DESeq(dds_wr, test = "LRT", reduced = ~1)


## Quick first looks
#What kind of comparisons did your chosen test conduct?
resultsNames(dds_wr)

# Total number of raw counts per sample
colSums(counts(dds_wr))

# Total number of normalized counts per sample
colSums(counts(dds_wr, normalized=T))

# Plot dispersion estimates (use this page for explanation https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html)
plotDispEsts(dds_wr)


## Explore the log2fold change
res_wr <- results(dds_wr, alpha = 0.01) # automatically selects the first comparison for LFC value change in resultsNames(dds_wr)
summary(res_wr)
#res_test <- results(dds_wr, contrast = c("Category", "Post_Flowering", "Pre_Flowering"), alpha = 0.01) # for Wald test
#summary(res_test)

#Use coefficients to shrink (appropriate for LRT)
res_shrunk_wr <- lfcShrink(dds_wr, res = res_wr, type = "ashr") # Using this method for now
summary(res_shrunk_wr)
#Use contrast to shrink (appropriate for Wald Test)
#res_shrunk <- lfcShrink(dds_wr, contrast=c("Category", "Post_Flowering", "Pre_Flowering"), res = res_test, type = "ashr")
#summary(res_shrunk)

# Plot MA plots (blue = significant)
DESeq2::plotMA(res_wr, ylim=c(-10,10))
DESeq2::plotMA(res_shrunk_wr, ylim=c(-10,10))

# Make a simple volcano plot
plot(res_shrunk_wr$log2FoldChange, -log10(res_shrunk_wr$padj))

# A more sophisticated volcano plot
res_shrunk_wr_DF <- data.frame(res_shrunk_wr)
res_shrunk_wr_DF <- res_shrunk_wr_DF %>% 
  dplyr::mutate(threshold = padj < 0.01 & abs(log2FoldChange) >= 0.58) 


# Volcano plot
ggplot(res_shrunk_wr_DF) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold), size = 1) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_color_manual(labels = c("FALSE", "TRUE"), values=c("grey", "#F8766D")) +
  scale_x_continuous(limits = c(-20,20)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

## Look at the results output
# If using LRT, the p-value doesn't correspond to log2foldchange, so it's not associated with the actual test
res_wr
res_shrunk_wr

# Plot the most significantly deferentially expressed gene
plot_count_wr <- plotCounts(dds_wr, gene = which.min(res_wr$padj), 
                            intgroup = c("Time","Category"), returnData = TRUE)

ggplot(plot_count_wr,
       aes(x = Time, y = count, color = Category, group = Category)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +  
  scale_color_manual(labels = c("Pre_Flowering", "Post_Flowering"), values=c("seagreen3", "magenta2")) +
  ylab("Normalized Count")

# Which gene is this?
res_wr_view <- results(dds_wr, tidy=TRUE)
res_wr_view <- as_tibble(res_wr_view)
res_wr_view <- res_wr_view %>% arrange(padj)

head(res_wr_view, 20)

# Top 30 significantly expressed genes
Top30_sig_df <- res_wr_view[1:30, ]
Top30_sig_list <- Top30_sig_df$row # this will be used later for Go analysis

# Most expressed gene, generally?
res_wr_view_exp<- res_wr_view %>% arrange(desc(baseMean))

head(res_wr_view_exp, 20)

# Top 30 most expressed genes
Top30_exp_df <- res_wr_view_exp[1:30, ]
Top30_exp_list <- Top30_exp_df$row # this will be used later for GO analysis

# Plot these by specifying the gene
plot_gene_wr <- plotCounts(dds_wr, gene = "g11025", 
                           intgroup = c("Time","Category"), returnData = TRUE)

ggplot(plot_gene_wr,
       aes(x = Time, y = count, color = Category, group = Category)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +  
  scale_color_manual(labels = c("Pre_Flowering", "Post_Flowering"), values=c("seagreen3", "magenta2")) +
  ylab("Normalized Count")




######## Using DEGreport::degPatterns to identify clusters of genes w shared expression profiles ########

## Let's get the data ready
# Subset the LRT results to return genes with padj < 0.01
sig_res_wr <- res_shrunk_wr %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(log2FoldChange >= 0.58 | log2FoldChange <= -0.58) %>%  #This is where things can change!
  filter(padj < 0.01)

# Save this list
write.csv(sig_res_wr, file="df_sig_all_wr_Test3_2.csv")

# Arrange in order of significance, optional
sig_res_wr <- sig_res_wr %>%
  arrange(padj)

# Get number of significant genes
summary(sig_res_wr)
nrow(sig_res_wr)

# Remember extracting the rlog matrix? This is performed above as well
rld_wr <- rlog(dds_wr, blind=TRUE)
rld_wr_mat <- assay(rld_wr)
head(rld_wr_mat)

# Obtain rlog values for those significant genes
cluster_rlog_wr <- rld_wr_mat[sig_res_wr$gene, ]

## degPatterns: Hierarchical clustering approach based on pair-wise correlations
#Make first colum of metadata row names
meta_pr2 <- meta_pr[,-1]
rownames(meta_pr2) <- meta_pr[,1]

#Run degPatterns
clusters_wr <- degPatterns(cluster_rlog_wr, metadata = meta_pr2, time="Time", col="Category", minc = 0)
#clusters_wr <- degPatterns(cluster_rlog_wr, metadata = meta_wr2, time="Time", col="Category", nClusters=22)

# What type of data structure is the `clusters` output?
class(clusters_wr)

# Let's see what is stored in the `df` component
head(clusters_wr$df)

#Plot
degPlotCluster(clusters_wr$normalized, "Time", col="Category")
#degPlotCluster(clusters_wr$normalized, "Time", color = "blue")

# Save the cluster assignments
cluster_groups_wr <- clusters_wr$df
write.csv(cluster_groups_wr, file="df_cluster_groups_sig_wr_test3_2.csv")


######## Using Photosynthetic genes to annotate everything
# Load cleaned version
ortho_clean <- read.csv(file = 'photosynthetic_genes_orthologs_longer.csv')

# Keep relevant column
ortho_anno <- ortho_clean %>%
  dplyr::select(Label, Gene_family, Block_description, Pathway, Cislong_ortho)

View(ortho_anno)

# Extract the Cislong_ortho column from ortho_anno
photo_genes <- ortho_anno$Cislong_ortho

# Create a new column in res_shrunk_wr_DF indicating if the gene is present in ortho_anno
res_shrunk_wr_DF_anno <- res_shrunk_wr_DF %>%
  rownames_to_column(var = "gene_id") %>%   # Convert row names to a column for easier manipulation
  mutate(photo_genes_TF = gene_id %in% photo_genes) %>%   # Create the new column
  column_to_rownames(var = "gene_id")   # Convert the gene_id column back to row names

# Create a new dataframe with the concatenated gene labels for each gene
labels_df <- ortho_anno %>%
  group_by(Cislong_ortho) %>%
  summarize(Label_concat = paste(Label, collapse = ", ")) 

# Join this dataframe with the original res_shrunk_wr_DF_anno to add the new column
res_shrunk_wr_DF_anno <- res_shrunk_wr_DF_anno %>%
  rownames_to_column(var = "gene_id") %>%  # Convert rownames to a column for the join
  left_join(labels_df, by = c("gene_id" = "Cislong_ortho")) %>%  # Join on gene_id and Cislong_ortho
  column_to_rownames(var = "gene_id")  # Convert gene_id back to rownames

# Now for pathway labels
pathway_df <- ortho_anno %>%
  group_by(Cislong_ortho) %>%
  summarize(pathway = dplyr::first(Pathway))

# Join this dataframe with the original res_shrunk_wr_DF_anno to add the new column
res_shrunk_wr_DF_anno <- res_shrunk_wr_DF_anno %>%
  rownames_to_column(var = "gene_id") %>%  # Convert rownames to a column for the join
  left_join(pathway_df, by = c("gene_id" = "Cislong_ortho")) %>%  # Join on gene_id and Cislong_ortho
  column_to_rownames(var = "gene_id")  # Convert gene_id back to rownames

# View the updated dataframe
head(res_shrunk_wr_DF_anno)
View(res_shrunk_wr_DF_anno)

## Volcano plot
ggplot(res_shrunk_wr_DF_anno) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold), size = 0.8, alpha = 0.5) +
  geom_point(data = subset(res_shrunk_wr_DF_anno, photo_genes_TF == TRUE),
             aes(x = log2FoldChange, y = -log10(padj)),
             colour = "blue", size = 1, shape = 21, stroke = 0.5) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_color_manual(labels = c("FALSE", "TRUE"), values=c("grey", "#F8766D")) +
  scale_x_continuous(limits = c(-20,20)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

# Calculate the min and max for x (log2FoldChange) and y (-log10(padj)) in the subset of interest
subset_data <- subset(res_shrunk_wr_DF_anno, pathway == "Carboxylation")

min_value_x <- min(subset_data$log2FoldChange, na.rm = TRUE)
max_value_x <- max(subset_data$log2FoldChange, na.rm = TRUE)
min_value_y <- min(-log10(subset_data$padj), na.rm = TRUE)
max_value_y <- max(-log10(subset_data$padj), na.rm = TRUE)

# Ensure some padding around the points
padding_x <- 0.1 * (max_value_x - min_value_x)
padding_y <- 0.1 * (max_value_y - min_value_y)

## More Volcano plots
ggplot(res_shrunk_wr_DF_anno) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold), size = 0.8, alpha = 0.5) +
  geom_point(data = subset(res_shrunk_wr_DF_anno, pathway == "Carboxylation"),
             aes(x = log2FoldChange, y = -log10(padj)),
             colour = "blue", size = 1, shape = 21, stroke = 0.5) +
  geom_text_repel(data = subset(res_shrunk_wr_DF_anno, pathway == "Carboxylation"),
                  aes(x = log2FoldChange, y = -log10(padj), label = Label_concat),
                  size = 3,
                  color = "blue",
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  max.overlaps = 30) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_color_manual(labels = c("FALSE", "TRUE"), values=c("grey", "#F8766D")) +
  scale_x_continuous(limits = c(min_value_x - padding_x, max_value_x + padding_x)) +
  scale_y_continuous(limits = c(min_value_y - padding_y, max_value_y + padding_y)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

# Just PEPC
ppc_genes <- c("g7768", "g11025", "g9957", "g22069")

ggplot(res_shrunk_wr_DF_anno) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold), size = 0.8, alpha = 0.5) +
  geom_point(data = subset(res_shrunk_wr_DF_anno, rownames(res_shrunk_wr_DF_anno) %in% ppc_genes),
             aes(x = log2FoldChange, y = -log10(padj)),
             colour = "blue", size = 1.5, shape = 21, stroke = 0.5) +
  geom_text_repel(data = subset(res_shrunk_wr_DF_anno, rownames(res_shrunk_wr_DF_anno) %in% ppc_genes),
                  aes(x = log2FoldChange, y = -log10(padj), label = Label_concat),
                  size = 4, 
                  color = "blue",
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = "grey50") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_color_manual(labels = c("FALSE", "TRUE"), values=c("grey", "#F8766D")) +
  scale_x_continuous(limits = c(min_value_x - padding_x, max_value_x + padding_x)) +
  scale_y_continuous(limits = c(min_value_y - padding_y, max_value_y + padding_y)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

# Just the genes we care about:
important_genes <- c("g2548", "g2034", "g11025", "g17299", "g24696", "g25079", "g22207")

# Calculate the min and max for x (log2FoldChange) and y (-log10(padj)) in the subset of interest
subset_data <- subset(res_shrunk_wr_DF_anno, rownames(res_shrunk_wr_DF_anno) %in% important_genes)

min_value_x <- min(subset_data$log2FoldChange, na.rm = TRUE)
max_value_x <- max(subset_data$log2FoldChange, na.rm = TRUE)
min_value_y <- min(-log10(subset_data$padj), na.rm = TRUE)
max_value_y <- max(-log10(subset_data$padj), na.rm = TRUE)

# Ensure some padding around the points
padding_x <- 0.1 * (max_value_x - min_value_x)
padding_y <- 0.1 * (max_value_y - min_value_y)

# Volcano plot
ggplot(res_shrunk_wr_DF_anno) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold), size = 0.8, alpha = 0.5) +
  geom_point(data = subset(res_shrunk_wr_DF_anno, rownames(res_shrunk_wr_DF_anno) %in% important_genes),
             aes(x = log2FoldChange, y = -log10(padj)),
             colour = "blue", size = 1.5, shape = 21, stroke = 0.5) +
  geom_text_repel(data = subset(res_shrunk_wr_DF_anno, rownames(res_shrunk_wr_DF_anno) %in% important_genes),
                  aes(x = log2FoldChange, y = -log10(padj), label = Label_concat),
                  size = 4, 
                  color = "blue",
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = "grey50") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_color_manual(labels = c("FALSE", "TRUE"), values=c("grey", "#F8766D")) +
  scale_x_continuous(limits = c(min_value_x - padding_x, max_value_x + padding_x)) +
  scale_y_continuous(limits = c(min_value_y - padding_y, max_value_y + padding_y)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))


## Let's save the data associated with the photosynthetic genes
# Which are the photosynthetic genes?
photo_wr <- res_shrunk_wr_DF_anno %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(photo_genes_TF == TRUE)

write.csv(photo_wr, file="df_photo_wr.csv")

# Which of the photosynthetic genes are significantly DE?
sig_photo_wr <- res_shrunk_wr_DF_anno %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(photo_genes_TF == TRUE) %>%  # First filter by photo_genes_TF
#  filter(log2FoldChange >= 0.58 | log2FoldChange <= -0.58) %>%  # Then filter by log2FoldChange
  filter(padj < 0.01)  # Finally, filter by adjusted p-value

write.csv(sig_photo_wr, file="df_sig_photo_wr_test3_2.csv")


## Cluster analysis with the photosynthetic genes
# Arrange in order of significance, optional
sig_photo_wr <- sig_photo_wr %>%
  arrange(padj)

# Get number of significant genes
summary(sig_photo_wr)
nrow(sig_photo_wr)

# Extracting the rlog matrix is performed above as well so skip if done
rld_wr <- rlog(dds_wr, blind=TRUE)
rld_wr_mat <- assay(rld_wr)
head(rld_wr_mat)

# Obtain rlog values for those significant genes
cluster_rlog_photo_wr <- rld_wr_mat[sig_photo_wr$gene, ]

#Make first colum of metadata row names
meta_wr2 <- meta_wr[,-1]
rownames(meta_wr2) <- meta_wr[,1]

#Run degPatterns (Hierarchical clustering approach based on pair-wise correlations)
clusters_photo_wr <- degPatterns(cluster_rlog_photo_wr, metadata = meta_pr2, time="Time", col="Category", minc=0)
#clusters_wr <- degPatterns(cluster_rlog_wr, metadata = meta_wr2, time="Time", col="Category", nClusters=22)

# What type of data structure is the `clusters` output?
class(clusters_photo_wr)

# Let's see what is stored in the `df` component
head(clusters_photo_wr$df)

# Plot
degPlotCluster(clusters_photo_wr$normalized, "Time", col="Category")
#degPlotCluster(clusters_wr$normalized, "Time", color = "blue")

# Save the cluster assignments
cluster_groups_photo_wr <- clusters_photo_wr$df

View(cluster_groups_photo_wr)
write.csv(cluster_groups_photo_wr, file="df_cluster_groups_sig_photo_wr_test3_3.csv")

# Join this dataframe with the original res_shrunk_wr_DF_anno to add the new column
sig_photo_wr_anno <- sig_photo_wr %>%
  left_join(cluster_groups_photo_wr, by = c("gene" = "genes"))

## Let's create a bar graph by counting up the number of genes
sig_photo_wr_anno_summary <- sig_photo_wr_anno %>%
  group_by(cluster, pathway) %>%
  summarise(gene_count = n()) %>%
  ungroup()

# Plot
ggplot(sig_photo_wr_anno_summary, aes(x = as.factor(cluster), y = gene_count, fill = pathway)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cluster", y = "Number of Genes", fill = "Pathway") +
  theme_minimal()

# Now, with mean expression
sig_photo_wr_anno_summary <- sig_photo_wr_anno %>%
  group_by(cluster, pathway) %>%
  summarise(total_baseMean = sum(baseMean, na.rm = TRUE)) %>%
  ungroup()

# Plot
ggplot(sig_photo_wr_anno_summary, aes(x = as.factor(cluster), y = total_baseMean, fill = pathway)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cluster", y = "Total Mean Expression", fill = "Pathway") +
  theme_minimal()



######## Using IPR & GO terms for annotation and enrichment

## Load annotation spreadsheet
pep_ipr <- read.csv(file = 'cislon.pep_IPR.tsv', sep = '\t')

pep_ipr_filtered <- pep_ipr %>%
  filter(grepl("^g\\d+\\.t\\d+$", Gene)) %>%
  separate(Gene, into = c("Gene", "Transcript"), sep = "\\.")

View(pep_ipr_filtered)

# Clean up more
pep_func_ipr_go <- pep_ipr_filtered %>%
  filter(!(Function == "-" & IPRcode == "-" & GOterm1 == "-")) %>%
  filter(!Function == "consensus disorder prediction") 

View(pep_func_ipr_go)

## Replace "-" under IPR description with terms in under Function, if possible
pep_func_ipr_go$Description <- ifelse(pep_func_ipr_go$Description == "-", 
                                      pep_func_ipr_go$Function, 
                                      pep_func_ipr_go$Description)

View(pep_func_ipr_go)

# Keep only unique conbinations 
df_disc_unique <- pep_func_ipr_go %>% distinct(Gene, Description, .keep_all = TRUE)

View(df_disc_unique)

# Clean up more
df_disc_unique <- df_disc_unique %>%
  filter(!Function == "Coil")

View(df_disc_unique)

# Keep one row per gene
model_priority <- c("Hamap", "PIRSF", "PANTHER", "FunFam", "SMART", "SFLD", 
                    "Pfam", "ProSitePatterns", "ProSiteProfiles", "NCBIfam", 
                    "CDD", "SUPERFAMILY", "Gene3D", "PRINTS")

# Create a new column in df_disc_unique that gives the priority of the model
df_disc_unique <- df_disc_unique %>%
  mutate(Model_Priority = match(Model, model_priority))

# Arrange the dataframe by Gene, Model_Priority (ascending), and Pvalue (ascending)
df_disc_unique <- df_disc_unique %>%
  arrange(Gene, Model_Priority, Pvalue)

# For each gene, keep only the first row (the one with the highest priority and lowest Pvalue)
df_final_unique <- df_disc_unique %>%
  group_by(Gene) %>%
  slice(1) %>%
  ungroup()

# Remove the Model_Priority column if you no longer need it
df_final_unique <- df_final_unique %>%
  select(-Model_Priority)

View(df_final_unique)

# Ge this from the photosynthetic anno pipeline from above
# Join this dataframe with the original res_shrunk_wr_DF_anno to add the new column
res_shrunk_wr_DF_FULL <- res_shrunk_wr_DF_anno %>%
  rownames_to_column(var = "gene_id") %>%  # Convert rownames to a column for the join
  left_join(df_final_unique %>% select(Gene, IPRcode, Description, Model, ID, Function), by = c("gene_id" = "Gene")) %>%
  column_to_rownames(var = "gene_id") 

# Which of the genes are significantly DE? Now, I get the annotated version
sig_photo_wr_FULL <- res_shrunk_wr_DF_FULL %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(threshold == TRUE)  # This does the selecting for pvalue and log2fold change

write.csv(sig_photo_wr_FULL, file="df_sig_all_wr_ANNO_test3-2.csv")




######## GO Enrichment analyses

## Filter out only the row with GO_IDs
pep_go_only <- pep_func_ipr_go %>%
  filter(!(GOterm1 == "-"))

View(pep_go_only)

# Pivot longer
pep_go_only_long <- pep_go_only %>%
  pivot_longer(cols = starts_with("GOterm"),
               names_to = "GOterm_type",
               values_to = "GOterm_value",
               values_drop_na = TRUE) %>%
  filter(GOterm_value != "" & !is.na(GOterm_value)) %>%
  mutate(GOterm_value = str_remove(GOterm_value, "\\(.*\\)")) %>%
  dplyr::select(-GOterm_type) %>% # Remove the helper column
  filter(!is.na(GOterm_value))

View(pep_go_only_long)

# Assume df is your dataframe with columns gene_id and go_id
df_go <- pep_go_only_long %>%
  dplyr::select(Gene,
                GOterm_value) %>%
  rename(gene_id = Gene) %>%
  rename(go_id = GOterm_value)

View(df_go)

df_go_unique <- df_go %>% distinct(gene_id, go_id, .keep_all = TRUE)

# Extract unique GO IDs
unique_go_ids <- unique(df_go_unique$go_id)

# Retrieve GO terms
go_terms <- AnnotationDbi::select(GO.db, keys = unique_go_ids, columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
# Or, filter for biological processes (BP)
bp_terms <- subset(go_terms, ONTOLOGY == "BP")

# Merge terms with your dataframe
df_go_with_terms <- merge(df_go_unique, bp_terms, by.x = "go_id", by.y = "GOID", all.x = TRUE)

View(df_go_with_terms)

# Create a new dataframe with the concatenated GO terms for each gene
GO_df <- df_go_with_terms %>%
  filter(!is.na(TERM)) %>%
  group_by(gene_id) %>%
  summarize(Term = paste(TERM, collapse = ", ")) 

View(GO_df)

# Join this dataframe with res_shrunk_wr_DF_FULL
res_shrunk_wr_DF_FULL_GO <- res_shrunk_wr_DF_FULL %>%
  rownames_to_column(var = "gene_id") %>%  # Convert rownames to a column for the join
  left_join(GO_df %>% select(gene_id, Term), by = c("gene_id" = "gene_id")) %>%
  column_to_rownames(var = "gene_id") 

View(res_shrunk_wr_DF_FULL_GO)

# Which of the genes are significantly DE? Now, I get the annotated version
sig_photo_wr_FULL_GO <- res_shrunk_wr_DF_FULL_GO %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(threshold == TRUE)  # This does the selecting for pvalue and log2fold change

write.csv(sig_photo_wr_FULL_GO, file="df_sig_all_wr_ANNO_test3-2.csv")


# Attach terms with most significant or most expressed list from before
# But I should really be doing this will all the sig expressed genes, and then match the clusters
#df_with_terms_subset <- df_with_terms[df_with_terms$gene_id %in% Top30_sig_list, ]
#write.csv(df_with_terms_subset, file="Top30_with_terms.csv")

# Assume deg_list is your list of differentially expressed genes
deg_list <- sig_res_wr$gene

# Prepare TERM2GENE and TERM2NAME
term2gene <- df_go_with_terms[, c("go_id", "gene_id")]
term2name <- bp_terms[, c("GOID", "TERM")]
colnames(term2name) <- c("ID", "Description")

# Check for discrepancies
term2gene_ids <- unique(term2gene$go_id)
term2name_ids <- unique(term2name$ID)

# making sure term2gene and term2name overlaps
missing_ids <- setdiff(term2gene_ids, term2name_ids)
if (length(missing_ids) > 0) {
  cat("Missing GO IDs in TERM2NAME:\n")
  print(missing_ids)
  # Filter TERM2GENE to include only GO IDs present in TERM2NAME
  term2gene <- term2gene[term2gene$go_id %in% term2name$ID, ]
} else {
  cat("All GO IDs in TERM2GENE are present in TERM2NAME.\n")
}

all_genes_present <- all(deg_list %in% pep_go_only_long$Gene)

if(all_genes_present) {
  cat("All genes in deg_list are present in TERM2GENE.\n")
} else {
  cat("Some genes in deg_list are not present in TERM2GENE. Missing genes:\n")
  print(setdiff(deg_list, pep_go_only_long$Gene))
}

# pulling out only the valid GO terms, if "Missing GO IDs in TERM2NAME above"
#valid_term2gene <- term2gene[term2gene$go_id %in% term2name$ID, ]
deg_with_go <- deg_list[deg_list %in% unique(term2gene$gene_id)]
missing_genes <- setdiff(deg_list, deg_with_go)

print(deg_with_go)

if (length(missing_genes) > 0) {
  cat("Genes in DEG list without GO IDs in TERM2NAME:\n")
  print(missing_genes)
} else {
  cat("All DEG list genes have GO IDs in TERM2NAME.\n")
}

#Pull out background genes

background_genes_df <- res_shrunk_wr %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

background_genes <- background_genes_df$gene

# Run enrichment analysis using DEG list
ego <- enricher(gene         = deg_with_go,
                TERM2GENE    = term2gene,
                TERM2NAME    = term2name,
                pAdjustMethod = "none",
                pvalueCutoff = 0.2,  # Adjust p-value cutoff if needed
                qvalueCutoff = 0.7,
                universe     = background_genes)

# Visualize results
summary(ego)
dotplot(ego)
barplot(ego, showCategory=50)

head(ego)

write.csv(ego, file="df_wr_GO_BP_test3-2.csv")
write.csv(ego, file="df_Test1_Overlap_GO_BP.csv")


## GSEA analyses
all_res_wr <- res_shrunk_wr %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

ranked_genes_df <- all_res_wr %>%
  arrange(desc(log2FoldChange))

ranked_genes <- setNames(ranked_genes_df$log2FoldChange, ranked_genes_df$gene)

gsea_result <- GSEA(ranked_genes,
                    TERM2GENE = term2gene,
                    TERM2NAME = term2name,
                    pvalueCutoff = 0.1,  # Adjust p-value cutoff
                    pAdjustMethod = "none",
                    verbose = FALSE)

# Visualize GSEA results
dotplot(gsea_result)

head(gsea_result)

write.csv(gsea_result, file="df_wr_GSEA_BP_test3-2.csv")

# GSEA plot for a specific pathway
gseaplot(gsea_result, geneSetID = gsea_result$ID[1])




######## Analyzing expression levels (NOT Differential expression) ########
# create a vector of gene lengths 
myco_wr_raw <- mycounts_raw %>%
  dplyr::select(Length,
                L9B6PM,
                L9B2PM,
                L9B10PM,
                L9B10AM,
                L8K6AM,
                L8K2PM,
                L8K2AM,
                L8K10AM,
                L8A6PM,
                L8A6AM,
                L8A2PM,
                L8A10PM,
                L8A10AM,
                L5I6PM,
                L5I6AM,
                L5I2PM,
                L5I2AM,
                L5I10PM,
                L5I10AM,
                L4L6PM,
                L4L6AM,
                L4L2PM,
                L4L2AM,
                L4L10PM,
                L4L10AM,
                L3I6PM,
                L3I6AM,
                L3I2PM,
                L3I2AM,
                L3I10PM,
                L3I10AM,
                L2Z6PM,
                L2Z6AM,
                L2Z2PM,
                L2Z2AM,
                L2Z10PM,
                L2Z10AM,
                L10O6PM,
                L10O6AM,
                L10O2PM,
                L10O2AM,
                L10O10PM,
                L10O10AM)
View(myco_wr_raw)


######## Using Photosynthetic genes to annotate everything
ortho_long <- read.csv(file = 'photosynthetic_genes_orthologs_longer.csv')


######## Calculate TPM
# Get Gene length
geneLengths <- as.vector(subset(myco_wr_raw, select = c(Length)))

# Divide each gene count per sample by gene length, normalize by 1 million, normalize by total counts per sample 
x <- myco_wr_raw / geneLengths
myco_wr_tpm <- t( t(x) * 1e6 / colSums(x) )

View(myco_wr_tpm) # Does "Length have the same number?

## OR (same thing as above but easier to see what's going on)
#Calculate RPK
#gene_length <- as.numeric(as.vector(myco_wr_raw$Length))
#rpk <- myco_wr_raw / (gene_length / 1000)

#Calculate the scaling factor
#scaling_factor <- colSums(rpk) / 1e6

# Calculate TPM
#tpm <- sweep(rpk, 2, scaling_factor, FUN = "/")
#View(tpm)

# Delete the length column
myco_wr_tpm <- myco_wr_tpm[ ,-1]


#### Script to convert Time to Zeitgeber time
convert_to_zeitgeber <- function(time_str) {
  time_map <- c("10AM" = 4, "2PM" = 8, "6PM" = 12, "10PM" = 16, "2AM" = 20, "6AM" = 24)
  return(time_map[time_str])
}


#### Subset PPC genes
ppc_genes <- c("g7768", "g11025", "g9957", "g22069")

ortho_long_subset <- ortho_long[ortho_long$Cislong_ortho %in% ppc_genes, ]
View(ortho_long_subset)

# Change WD
setwd("~/Dropbox/Yale/Research/RNAseq/TPM_graphs_WR_NEW")

# Step 1: Extract count data for the specific gene
gene_id <- "g11025"  # Replace with your gene of interest
gene_df <- myco_wr_tpm[rownames(myco_wr_tpm) == gene_id, , drop = FALSE]

# Ensure the gene exists in the dataset
if (nrow(gene_df) == 0) {
  stop(paste("Gene ID", gene_id, "not found in myco_wr_tpm"))
}

# Step 2: Reshape the data and change the first column name
gene_df <- as.data.frame(t(gene_df))
gene_df$SampleID <- rownames(gene_df)

colnames(gene_df)[1] = "counts"

# Step 3: Merge the count data with the metadata-
# Now,`merged_data` has the `count`, `Time`, and `Category` columns along with other metadata
merged_data <- merge(gene_df, meta_wr, by = "SampleID")

# Step 4: Convert Time to Zeitgeber_Time
merged_data$Zeitgeber_Time <- sapply(merged_data$Time, convert_to_zeitgeber)

# Step 5: Make Summary of the data to compute means and SD
summary_df <- merged_data %>%
  group_by(Zeitgeber_Time, Category) %>%
  summarise(mean_counts = mean(counts), se_counts = sd(counts)/sqrt(n()))

# Step 6: Plotting
p <- ggplot(summary_df, aes(x = Zeitgeber_Time, y = mean_counts, color = Category, group = Category)) + 
  geom_errorbar(aes(ymin = mean_counts - se_counts, ymax = mean_counts + se_counts), width = 0.5, linewidth = 2) +
  geom_line(linewidth = 2) +  
  scale_color_manual(labels = c("Pre-Flowering", "Post-Flowering"), values=c("seagreen3", "magenta2")) +
  ylab("TPM") +
  xlab("Zeitgeber Time (hours)") +
  scale_x_continuous(breaks = seq(4, 24, by = 4), limits = c(3, 24)) +
  annotate("rect", xmin = 14, xmax = 24, ymin = -Inf, ymax =Inf, fill = "lightgrey", alpha = 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size=20)) +
  ggtitle(paste0("PPC-1E1c | ", gene_id))

p

p <- ggplot(merged_data, aes(x = Zeitgeber_Time, y = counts, color = Category, group = Category)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_color_manual(labels = c("Pre_Flowering", "Post_Flowering"), values=c("seagreen3", "magenta2")) +
  ylab("TPM") +
  xlab("Zeitgeber Time (hours)") +
  scale_x_continuous(breaks = seq(4, 24, by = 4), limits = c(3, 25)) +
  annotate("rect", xmin = 14, xmax = 24, ymin = -Inf, ymax = Inf, fill = "lightgrey", alpha = 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_blank()) +
  ggtitle(paste0("Gene:", gene_id))



#### Loop through each gene in the subset or the full dataset

# Change WD
setwd("~/Dropbox/Yale/Research/RNAseq/TPM_graphs_WR_axis")

# Loop
for (i in 1:nrow(ortho_long)) {
  gene_id <- ortho_long$Cislong_ortho[i]
  label <- ortho_long$Label[i]
  
  # Step 1: Extract count data for the specific gene
  gene_df <- myco_wr_tpm[rownames(myco_wr_tpm) == gene_id, , drop = FALSE]
  
  # Ensure the gene exists in the dataset
  if (nrow(gene_df) == 0) {
    warning(paste("Gene ID", gene_id, "not found in myco_wr_tpm Skipping..."))
    next
  }
  
  # Step 2: Reshape the data using pivot_longer
  gene_df <- as.data.frame(t(gene_df))
  gene_df$SampleID <- rownames(gene_df)
  
  colnames(gene_df)[1] = "counts"
  
  # Step 3: Merge the count data with the metadata
  merged_data <- merge(gene_df, meta_wr, by = "SampleID")
  
  # Convert Time to Zeitgeber time
  merged_data$Zeitgeber_Time <- sapply(merged_data$Time, convert_to_zeitgeber)
  
  # Summarize data by Zeitgeber time and Category
  summary_df <- merged_data %>%
    group_by(Zeitgeber_Time, Category) %>%
    summarise(mean_counts = mean(counts), se_counts = sd(counts)/sqrt(n()))
  
  # Step 4: Plotting
  p <- ggplot(summary_df, aes(x = Zeitgeber_Time, y = mean_counts, color = Category, group = Category)) + 
    geom_errorbar(aes(ymin = mean_counts - se_counts, ymax = mean_counts + se_counts), width = 0.5, linewidth = 2.5) +
    geom_line(linewidth = 2.5) +  
    scale_color_manual(labels = c("Pre_Flowering", "Post_Flowering"), values=c("seagreen3", "magenta2")) + #43CD80, #EE00EE
    ylab("TPM") +
    xlab("Zeitgeber time (hours)") +
    scale_x_continuous(breaks = seq(4, 24, by = 4), limits = c(3, 24)) +
    annotate("rect", xmin = 14, xmax = 24, ymin = -Inf, ymax = Inf, fill = "lightgrey", alpha = 0.3) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size=22)) +
    ggtitle(paste0(label, " | ", gene_id))
  
  # Save the plot
  ggsave(filename = paste0(label, "_", gene_id, ".png"), plot = p)
}



