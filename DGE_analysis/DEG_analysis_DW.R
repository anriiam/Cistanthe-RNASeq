rm(list=ls())

# Commented out because I've already installed the program
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#  BiocManager::install(version = "3.19")
#BiocManager::install("DESeq2")
#BiocManager::install("maSigPro")
#BiocManager::install("DEGreport")


################ RNA DGE analysis Pipeline ################
# Anri Chomentowska August 2024


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
library(WGCNA)
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


#### Subset dataframes into the two experiments
## Subset (1) Drought v Watered experiment
meta_dw <- metadata %>%
  dplyr::filter(Daylength == "Short")

myco_dw <- mycounts %>%
  dplyr::select(S9A6PM,
                S9A6AM,
                S9A2AM,
                S9A10PM,
                S9A10AM,
                S6H6PM,
                S6H6AM,
                S6H2PM,
                S6H2AM,
                S6H10PM,
                S6H10AM,
                S4J6PM,
                S4J6AM,
                S4J2PM,
                S4J2AM,
                S4J10PM,
                S4J10AM,
                S4H6PM,
                S4H2PM,
                S4H10PM,
                S4H10AM,
                S4B6PM,
                S4B6AM,
                S4B2PM,
                S4B2AM,
                S4B10PM,
                S4B10AM,
                S12N6PM,
                S12N6AM,
                S12N2PM,
                S12N2AM,
                S12N10PM,
                S12N10AM,
                S11L6PM,
                S11L6AM,
                S11L2PM,
                S11L2AM,
                S11L10PM,
                S11L10AM,
                S10F6PM,
                S10F6AM,
                S10F2PM,
                S10F2AM,
                S10F10AM)

View(meta_dw)
View(myco_dw)

# Make sure the IDs are identical
names(myco_dw)
meta_dw$SampleID

names(myco_dw)==meta_dw$SampleID # Should say TRUE for all
all(names(myco_dw)==meta_dw$SampleID) # Should say TRUE once

# Set correct order for time course and how data will be shown
meta_dw$Time <- factor(meta_dw$Time, levels = c("10AM","2PM","6PM","10PM","2AM","6AM"))
meta_dw$Treatment <- factor(meta_dw$Treatment, levels = c("Watered", "Drought"))



######## QC & DESeq2 analyses ########

######## Experiment 1: Drought vs watered comparison
# Let's find significantly expressed genes with respect to time and treatment

#### Create DESeq2Dataset object
# Ignore the warning that comes up when it says "converting counts to integer mode..."
dds_dw <- DESeqDataSetFromMatrix(countData = round(myco_dw), 
                                 colData = meta_dw, 
                                 design = ~ Time + Treatment, 
                                 tidy = FALSE)
dds_dw # Take a look at the DESeq2 data structure

#### QC steps
# Perform the median of ratios method of normalization
# (This is also automatically performed by the DESeq() function, but you can take this output for use with other packages)
dds_dw <- estimateSizeFactors(dds_dw) # calculate normalization factor per sample (if significantly off from 1, then this means)

sizeFactors(dds_dw) #if very off from 1, then this means extreme outliers may exist

# We could tell that sample S6H2AM is an extreme outlier, looking along with the MultiQC output. Let's filter that out and redo the analysis
myco_dw <- myco_dw %>% dplyr::select(-S6H2AM)
meta_dw <- meta_dw %>% filter(SampleID != "S6H2AM")

View(myco_dw)
View(meta_dw)

# Make sure the IDs are identical, again
names(myco_dw)
meta_dw$SampleID

names(myco_dw)==meta_dw$SampleID # Should say TRUE for all
all(names(myco_dw)==meta_dw$SampleID) # Should say TRUE once

# Now rerun the DESeqMatrix and size factor calculation steps
dds_dw <- DESeqDataSetFromMatrix(countData = round(myco_dw), 
                                 colData = meta_dw, 
                                 design = ~ Time + Treatment, 
                                 tidy = FALSE)
dds_dw 
dds_dw <- estimateSizeFactors(dds_dw) # calculate normalization factor per sample (if significantly off from 1, then this means)
sizeFactors(dds_dw)


# retrieve the normalized counts matrix
normalized_counts_dw <- counts(dds_dw, normalized=TRUE)
write.table(normalized_counts_dw, file="normalized_counts_DW_NEW.txt", sep="\t", quote=F, col.names=NA) # Use this for analyses with other packages!
write.csv(meta_dw, file="metadata_DW_NEW.csv")

## Regularized log transformation (rlog) and Variance stabilizing transformation (vst) for data visualization
rld_dw <- rlog(dds_dw, blind=TRUE)
vst_dw <- vst(dds_dw, blind=TRUE)

# Extract the rlog and vst matrix from the object
rld_dw_mat <- assay(rld_dw)
vst_dw_mat <- assay(vst_dw)

# Compute pairwise correlation values
rld_dw_cor <- cor(rld_dw_mat)    ## cor() is a base R function
head(rld_dw_cor)   ## check the output of cor(), make note of the rownames and colnames

vst_dw_cor <- cor(vst_dw_mat)    ## cor() is a base R function
head(vst_dw_cor)   ## check the output of cor(), make note of the rownames and colnames

# Make first colum of metadata row names
meta_dw2 <- meta_dw[,-1]
rownames(meta_dw2) <- meta_dw[,1]
head(meta_dw2)

## Plot heatmap
pheatmap(rld_dw_cor, col = hmcol)
pheatmap(vst_dw_cor)

# Or, select the top 20 genes with the highest mean expression for heatmap
select_dw <- order(rowMeans(counts(dds_dw,normalized=TRUE)), decreasing=TRUE)[1:20]
anno_dw <- as.data.frame(colData(dds_dw)[,c("Time","Treatment")])

anno_dw$Time <- factor(anno_dw$Time, levels = c("10AM","2PM","6PM","10PM","2AM","6AM"))
anno_dw$Treatment <- factor(anno_dw$Treatment, levels = c("Watered", "Drought"))

# Order the columns in the matrix based on Treatment and then Time
ordered_cols <- order(anno_dw$Treatment, anno_dw$Time)
anno_dw_ordered <- anno_dw[ordered_cols, ]

rld_dw_mat_ordered <- rld_dw_mat[, ordered_cols]
vst_dw_mat_ordered <- vst_dw_mat[, ordered_cols]

# Plot the heatmap of for the selected top 20 genes to visualize
pheatmap(rld_dw_mat_ordered[select_dw,], cluster_rows=FALSE, col = hmcol,
         cluster_cols=FALSE, annotation_col=anno_dw_ordered)

pheatmap(vst_dw_mat_ordered[select_dw,], cluster_rows=FALSE, col = hmcol,
         cluster_cols=FALSE, annotation_col=anno_dw_ordered)

## Now, Plot a PCA 
plotPCA(rld_dw, intgroup="Treatment")
plotPCA(vst_dw, intgroup="Treatment")

# Different ways to explore more about the PCs (https://tavareshugo.github.io/data-carpentry-rnaseq/03_rnaseq_pca.html)
pca_dw <- prcomp(t(rld_dw_mat))
#pca_dw <- prcomp(t(vst_dw_mat)) #if you want to use vst
View(pca_dw)
df_dw <- cbind(meta_dw, pca_dw$x) # You can explore many ways to do this

pc_loadings_dw <- pca_dw$rotation
pc_loadings_dw <- pc_loadings_dw %>% 
  as_tibble(rownames = "gene")
View(pc_loadings_dw)

top_genes_dw <- pc_loadings_dw %>% 
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

top_genes_dw #list of the top genes (based on PC loadings)

# Filter the top loading genes' PC matrix and visualize the direction and magnitude of the eigenvalues
top_loadings_dw <- pc_loadings_dw %>% 
  filter(gene %in% top_genes_dw)

ggplot(top_loadings_dw) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "blue") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))

#another way to visualize and annotate PCA
pca_data_dw <- plotPCA(rld_dw, intgroup=c("Treatment","Time"), returnData=TRUE)
View(pca_data_dw)

percentVar_dw <- 100*round(attr(pca_data_dw, "percentVar"),2)               
makeLab <- function(x,pc) paste0("PC",pc,": ",x,"% variance")

ggplot(pca_data_dw) + geom_point(aes(x=PC1, y=PC2, color = Treatment, shape = Time), size=3) +
  xlab(makeLab(percentVar_dw[1],1)) + ylab(makeLab(percentVar_dw[2],2)) + 
  scale_color_manual(labels = c("Watered", "Drought"), values=c("#619CFF", "#F8766D")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_blank()) +
  ggtitle("PCA (rlog)")

#vst
pca_data_dw <- plotPCA(vst_dw, intgroup=c("Treatment","Time"), returnData=TRUE)
View(pca_data_dw)

percentVar_dw <- 100*round(attr(pca_data_dw, "percentVar"),2)               
makeLab <- function(x,pc) paste0("PC",pc,": ",x,"% variance")

ggplot(pca_data_dw) + geom_point(aes(x=PC1, y=PC2, color = Treatment, shape = Time), size=3) +
  xlab(makeLab(percentVar_dw[1],1)) + ylab(makeLab(percentVar_dw[2],2)) + 
  scale_color_manual(labels = c("Watered", "Drought"), values=c("#619CFF", "#F8766D")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_blank()) +
  ggtitle("PCA (vst)")


######## Differential gene expression (DEG) steps

## TEST 1 (testing just Treatment effects via LRT, for any differences over multiple time points)
dds_dw <- DESeqDataSetFromMatrix(countData = round(myco_dw), 
                                 colData = meta_dw, 
                                 design = ~ Time + Treatment, 
                                 tidy=FALSE)

dds_dw <- DESeq(dds_dw, test = "LRT", reduced = ~ Time)

# TEST 2 (p-value calculated from difference in treatment: LRT or Wald Test
dds_dw <- DESeqDataSetFromMatrix(countData = round(myco_dw),
                                 colData = meta_dw, 
                                 design = ~ Treatment, 
                                 tidy=FALSE)

dds_dw <- DESeq(dds_dw, test = "LRT", reduced = ~1) #LRT
dds_dw <- DESeq(dds_dw) #Wald test

dds_dw

# TEST 3.1 (p-value calculated from effect of Time, but between each treatment): Drought
meta_dr <- meta_dw %>%
  dplyr::filter(Treatment == "Drought")

myco_dr <- myco_dw %>%
  dplyr::select(S9A6PM,
                S9A6AM,
                S9A2AM,
                S9A10PM,
                S9A10AM,
                S4H6PM,
                S4H2PM,
                S4H10PM,
                S4H10AM,
                S4B6PM,
                S4B6AM,
                S4B2PM,
                S4B2AM,
                S4B10PM,
                S4B10AM,
                S11L6PM,
                S11L6AM,
                S11L2PM,
                S11L2AM,
                S11L10PM,
                S11L10AM)

dds_dw <- DESeqDataSetFromMatrix(countData = round(myco_dr), 
                                 colData = meta_dr, 
                                 design = ~ Time, 
                                 tidy=FALSE)
dds_dw <- DESeq(dds_dw, test = "LRT", reduced = ~1)

# TEST 3.2 (p-value calculated from effect of Time, but between each treatment): Watered
meta_wa <- meta_dw %>%
  dplyr::filter(Treatment == "Watered")

myco_wa <- myco_dw %>%
  dplyr::select(S6H6PM,
                S6H6AM,
                S6H2PM,
                S6H10PM,
                S6H10AM,
                S4J6PM,
                S4J6AM,
                S4J2PM,
                S4J2AM,
                S4J10PM,
                S4J10AM,
                S12N6PM,
                S12N6AM,
                S12N2PM,
                S12N2AM,
                S12N10PM,
                S12N10AM,
                S10F6PM,
                S10F6AM,
                S10F2PM,
                S10F2AM,
                S10F10AM)

dds_dw <- DESeqDataSetFromMatrix(countData = round(myco_wa), 
                                 colData = meta_wa, 
                                 design = ~ Time, 
                                 tidy=FALSE)
dds_dw <- DESeq(dds_dw, test = "LRT", reduced = ~1)


## Quick first looks
#What kind of comparisons did your chosen test conduct?
resultsNames(dds_dw)

# Total number of raw counts per sample
colSums(counts(dds_dw))

# Total number of normalized counts per sample
colSums(counts(dds_dw, normalized=T))

# Plot dispersion estimates (use this page for explanation https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html)
plotDispEsts(dds_dw)


## Explore the log2fold change
res_dw <- results(dds_dw, alpha = 0.01) # automatically selects the first comparison for LFC value change in resultsNames(dds_dw)
#res_test <- results(dds_dw, contrast = c("Treatment", "Drought", "Watered"), alpha = 0.01) # for Wald test

#Use coefficients to shrink (appropriate for LRT)
res_shrunk_dw <- lfcShrink(dds_dw, res = res_dw, type = "ashr") # Using this method for now
#Use contrast to shrink (appropriate for Wald Test)
#res_shrunk <- lfcShrink(dds_dw, contrast=c("Treatment", "Drought", "Watered"), res = res_test, type = "ashr")

# Plot MA plots (blue = significant)
DESeq2::plotMA(res_dw, ylim=c(-10,10))
DESeq2::plotMA(res_shrunk_dw, ylim=c(-10,10))

# None of the significant genes should change between the two
summary(res_dw)
summary(res_shrunk_dw)
#summary(res_test)
#summary(res_shrunk)

# Make a simple volcano plot
plot(res_shrunk_dw$log2FoldChange, -log10(res_shrunk_dw$padj))

# A more sophisticated volcano plot
res_shrunk_dw_DF <- data.frame(res_shrunk_dw)
res_shrunk_dw_DF <- res_shrunk_dw_DF %>% 
  dplyr::mutate(threshold = padj < 0.01 & abs(log2FoldChange) >= 0.58) 


# Volcano plot
ggplot(res_shrunk_dw_DF) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold), size = 1, alpha = 0.5) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_color_manual(labels = c("FALSE", "TRUE"), values=c("grey", "#F8766D")) +
  #scale_y_continuous(limits = c(0,50)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

## Look at the results output
# If using LRT, the p-value doesn't correspond to log2foldchange, so it's not associated with the actual test
res_dw
res_shrunk_dw

# Plot the most significantly deferentially expressed gene
plot_count_dw <- plotCounts(dds_dw, gene = which.min(res_dw$padj), 
                            intgroup = c("Time","Treatment"), returnData = TRUE)

ggplot(plot_count_dw,
       aes(x = Time, y = count, color = Treatment, group = Treatment)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +  
  scale_color_manual(labels = c(#"Watered",
    "Drought"), values=c(#"#619CFF",
    "#F8766D")) +
  ylab("Normalized Count")

# Which gene is this?
res_dw_view <- results(dds_dw, tidy=TRUE)
res_dw_view <- as_tibble(res_dw_view)
res_dw_view <- res_dw_view %>% arrange(padj)

head(res_dw_view, 20)

# Top 30 significantly expressed genes
Top30_sig_df <- res_dw_view[1:30, ]
Top30_sig_list <- Top30_sig_df$row # this will be used later for Go analysis

# Most expressed gene, generally?
res_dw_view_exp<- res_dw_view %>% arrange(desc(baseMean))

head(res_dw_view_exp, 20)

# Top 30 most expressed genes
Top30_exp_df <- res_dw_view_exp[1:30, ]
Top30_exp_list <- Top30_exp_df$row # this will be used later for GO analysis

# Plot these by specifying the gene, make sure it's the same plot from above
plot_gene_dw <- plotCounts(dds_dw, gene = "g11025", 
                    intgroup = c("Time","Treatment"), returnData = TRUE)

ggplot(plot_gene_dw,
       aes(x = Time, y = count, color = Treatment, group = Treatment)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +  
  scale_color_manual(labels = c("Watered","Drought"), values=c("#619CFF","#F8766D")) +
  ylab("Normalized Count")




######## Using DEGreport::degPatterns to identify clusters of genes w shared expression profiles ########

## Let's get the data ready
# Subset the LRT results to return genes with padj < 0.01
sig_res_dw <- res_shrunk_dw %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(log2FoldChange >= 0.58 | log2FoldChange <= -0.58) %>%  #This is where things can change!
  filter(padj < 0.01)

# Save this list
write.csv(sig_res_dw, file="df_sig_all_dw_Test2New.csv")

# Arrange in order of significance, optional
sig_res_dw <- sig_res_dw %>%
  arrange(padj)

# Get number of significant genes
summary(sig_res_dw)
nrow(sig_res_dw)

# Remember extracting the rlog matrix? This is performed above as well
rld_dw <- rlog(dds_dw, blind=TRUE)
rld_dw_mat <- assay(rld_dw)
head(rld_dw_mat)

# Obtain rlog values for those significant genes
cluster_rlog_dw <- rld_dw_mat[sig_res_dw$gene, ]

## degPatterns: Hierarchical clustering approach based on pair-wise correlations
#Make first colum of metadata row names
meta_wa2 <- meta_wa[,-1]
rownames(meta_wa2) <- meta_wa[,1]

#Run degPatterns
clusters_dw <- degPatterns(cluster_rlog_dw, metadata = meta_wa2, time="Time", col="Treatment", minc = 0)
#clusters_dw <- degPatterns(cluster_rlog_dw, metadata = meta_dw2, time="Time", col="Treatment", nClusters=22, minc > 0)

# What type of data structure is the `clusters` output?
class(clusters_dw)

# Let's see what is stored in the `df` component
head(clusters_dw$df)

#Plot
degPlotCluster(clusters_dw$normalized, "Time", col="Treatment")

# Save the cluster assignments
cluster_groups_dw <- clusters_dw$df
write.csv(cluster_groups_dw, file="df_cluster_groups_sig_all_dw_test3-2.csv")


######## Using Photosynthetic genes to annotate everything
# Load cleaned version
ortho_clean <- read.csv(file = 'photosynthetic_genes_orthologs_longer.csv')

# Keep relevant column
ortho_anno <- ortho_clean %>%
  dplyr::select(Label, Gene_family, Block_description, Pathway, Cislong_ortho)

View(ortho_anno)

# Extract the Cislong_ortho column from ortho_anno
photo_genes <- ortho_anno$Cislong_ortho

# Also done above
res_shrunk_dw_DF <- data.frame(res_shrunk_dw)
res_shrunk_dw_DF <- res_shrunk_dw_DF %>% 
  dplyr::mutate(threshold = padj < 0.01 & abs(log2FoldChange) >= 0.58) 

# Create a new column in res_shrunk_dw_DF indicating if the gene is present in ortho_anno
res_shrunk_dw_DF_anno <- res_shrunk_dw_DF %>%
  rownames_to_column(var = "gene_id") %>%   # Convert row names to a column for easier manipulation
  mutate(photo_genes_TF = gene_id %in% photo_genes) %>%   # Create the new column
  column_to_rownames(var = "gene_id")   # Convert the gene_id column back to row names

# Create a new dataframe with the concatenated gene labels for each gene
labels_df <- ortho_anno %>%
  group_by(Cislong_ortho) %>%
  summarize(Label_concat = paste(Label, collapse = ", ")) 

# Join this dataframe with the original res_shrunk_dw_DF_anno to add the new column
res_shrunk_dw_DF_anno <- res_shrunk_dw_DF_anno %>%
  rownames_to_column(var = "gene_id") %>%  # Convert rownames to a column for the join
  left_join(labels_df, by = c("gene_id" = "Cislong_ortho")) %>%  # Join on gene_id and Cislong_ortho
  column_to_rownames(var = "gene_id")  # Convert gene_id back to rownames

# Now for pathway labels
pathway_df <- ortho_anno %>%
   group_by(Cislong_ortho) %>%
   summarize(pathway = dplyr::first(Pathway))

# Join this dataframe with the original res_shrunk_dw_DF_anno to add the new column
res_shrunk_dw_DF_anno <- res_shrunk_dw_DF_anno %>%
  rownames_to_column(var = "gene_id") %>%  # Convert rownames to a column for the join
  left_join(pathway_df, by = c("gene_id" = "Cislong_ortho")) %>%  # Join on gene_id and Cislong_ortho
  column_to_rownames(var = "gene_id")  # Convert gene_id back to rownames

# View the updated dataframe
head(res_shrunk_dw_DF_anno)
View(res_shrunk_dw_DF_anno)

## Volcano plot
ggplot(res_shrunk_dw_DF_anno) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold), size = 0.8, alpha = 0.5) +
  geom_point(data = subset(res_shrunk_dw_DF_anno, photo_genes_TF == TRUE),
             aes(x = log2FoldChange, y = -log10(padj)),
             colour = "blue", size = 1, shape = 21, stroke = 0.5) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_color_manual(labels = c("FALSE", "TRUE"), values=c("grey", "#F8766D")) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

# Calculate the min and max for x (log2FoldChange) and y (-log10(padj)) in the subset of interest
subset_data <- subset(res_shrunk_dw_DF_anno, pathway == "Carboxylation")

min_value_x <- min(subset_data$log2FoldChange, na.rm = TRUE)
max_value_x <- max(subset_data$log2FoldChange, na.rm = TRUE)
min_value_y <- min(-log10(subset_data$padj), na.rm = TRUE)
max_value_y <- max(-log10(subset_data$padj), na.rm = TRUE)

# Ensure some padding around the points
padding_x <- 0.1 * (max_value_x - min_value_x)
padding_y <- 0.1 * (max_value_y - min_value_y)

## More Volcano plots
ggplot(res_shrunk_dw_DF_anno) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold), size = 0.8, alpha = 0.5) +
  geom_point(data = subset(res_shrunk_dw_DF_anno, pathway == "Carboxylation"),
             aes(x = log2FoldChange, y = -log10(padj)),
             colour = "blue", size = 1, shape = 21, stroke = 0.5) +
  geom_text_repel(data = subset(res_shrunk_dw_DF_anno, pathway == "Carboxylation"),
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
#  scale_x_continuous(limits = c(min_value_x - padding_x, max_value_x + padding_x)) +
#  scale_y_continuous(limits = c(min_value_y - padding_y, max_value_y + padding_y)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

# Just PEPC
ppc_genes <- c("g7768", "g11025", "g9957", "g22069")

ggplot(res_shrunk_dw_DF_anno) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold), size = 0.8, alpha = 0.5) +
  geom_point(data = subset(res_shrunk_dw_DF_anno, rownames(res_shrunk_dw_DF_anno) %in% ppc_genes),
             aes(x = log2FoldChange, y = -log10(padj)),
             colour = "blue", size = 1.5, shape = 21, stroke = 0.5) +
  geom_text_repel(data = subset(res_shrunk_dw_DF_anno, rownames(res_shrunk_dw_DF_anno) %in% ppc_genes),
                  aes(x = log2FoldChange, y = -log10(padj), label = Label_concat),
                  size = 4, 
                  color = "blue",
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = "grey50") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_color_manual(labels = c("FALSE", "TRUE"), values=c("grey", "#F8766D")) +
#  scale_x_continuous(limits = c(min_value_x - padding_x, max_value_x + padding_x)) +
#  scale_y_continuous(limits = c(min_value_y - padding_y, max_value_y + padding_y)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

# Just the genes we care about:
important_genes <- c("g2548", "g2034", "g11025", "g17299", "g24696", "g25079", "g22207")

# Calculate the min and max for x (log2FoldChange) and y (-log10(padj)) in the subset of interest
subset_data <- subset(res_shrunk_dw_DF_anno, rownames(res_shrunk_dw_DF_anno) %in% important_genes)

min_value_x <- min(subset_data$log2FoldChange, na.rm = TRUE)
max_value_x <- max(subset_data$log2FoldChange, na.rm = TRUE)
min_value_y <- min(-log10(subset_data$padj), na.rm = TRUE)
max_value_y <- max(-log10(subset_data$padj), na.rm = TRUE)

# Ensure some padding around the points
padding_x <- 0.1 * (max_value_x - min_value_x)
padding_y <- 0.1 * (max_value_y - min_value_y)

# Volcano plot
ggplot(res_shrunk_dw_DF_anno) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold), size = 0.8, alpha = 0.5) +
  geom_point(data = subset(res_shrunk_dw_DF_anno, rownames(res_shrunk_dw_DF_anno) %in% important_genes),
             aes(x = log2FoldChange, y = -log10(padj)),
             colour = "blue", size = 1.5, shape = 21, stroke = 0.5) +
  geom_text_repel(data = subset(res_shrunk_dw_DF_anno, rownames(res_shrunk_dw_DF_anno) %in% important_genes),
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
photo_dw <- res_shrunk_dw_DF_anno %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(photo_genes_TF == TRUE)

#write.csv(photo_dw, file="df_photo_dw.csv")

# Which of the photosynthetic genes are significantly DE?
sig_photo_dw <- res_shrunk_dw_DF_anno %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(photo_genes_TF == TRUE) %>%  # First filter by photo_genes_TF
#  filter(log2FoldChange >= 0.58 | log2FoldChange <= -0.58) %>%  # Then filter by log2FoldChange
  filter(padj < 0.01)  # Finally, filter by adjusted p-value

write.csv(sig_photo_dw, file="df_sig_photo_dw_test3-2.csv")


## Cluster analysis with the photosynthetic genes
# Arrange in order of significance, optional
sig_photo_dw <- sig_photo_dw %>%
  arrange(padj)

# Get number of significant genes
summary(sig_photo_dw)
nrow(sig_photo_dw)

# Remember extracting the rlog matrix? This is performed above as well so skip if done
rld_dw <- rlog(dds_dw, blind=TRUE)
rld_dw_mat <- assay(rld_dw)
head(rld_dw_mat)

# Obtain rlog values for those significant genes
cluster_rlog_photo_dw <- rld_dw_mat[sig_photo_dw$gene, ]

#Make first colum of metadata row names, also done above
meta_dr2 <- meta_dr[,-1]
rownames(meta_dr2) <- meta_dr[,1]

#Run degPatterns (Hierarchical clustering approach based on pair-wise correlations)
clusters_photo_dw <- degPatterns(cluster_rlog_photo_dw, metadata = meta_wa2, time="Time", col="Treatment", minc=0)
#clusters_dw <- degPatterns(cluster_rlog_dw, metadata = meta_dw2, time="Time", col="Treatment", nClusters=22)

# What type of data structure is the `clusters` output?
class(clusters_photo_dw)

# Let's see what is stored in the `df` component
head(clusters_photo_dw$df)

# Plot
degPlotCluster(clusters_photo_dw$normalized, "Time", col="Treatment")
#degPlotCluster(clusters_dw$normalized, "Time", color = "blue")

# Save the cluster assignments
cluster_groups_photo_dw <- clusters_photo_dw$df

View(cluster_groups_photo_dw)
write.csv(cluster_groups_photo_dw, file="df_cluster_groups_sig_photo_dw_test3-2.csv")

# Join this dataframe with the original res_shrunk_dw_DF_anno to add the new column
sig_photo_dw_anno <- sig_photo_dw %>%
  left_join(cluster_groups_photo_dw, by = c("gene" = "genes"))

## Let's create a bar graph by counting up the number of genes
sig_photo_dw_anno_summary <- sig_photo_dw_anno %>%
  group_by(cluster, pathway) %>%
  summarise(gene_count = n()) %>%
  ungroup()

# Plot
ggplot(sig_photo_dw_anno_summary, aes(x = as.factor(cluster), y = gene_count, fill = pathway)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cluster", y = "Number of Genes", fill = "Pathway") +
  theme_minimal()

# Now, with mean expression
sig_photo_dw_anno_summary <- sig_photo_dw_anno %>%
  group_by(cluster, pathway) %>%
  summarise(total_baseMean = sum(baseMean, na.rm = TRUE)) %>%
  ungroup()

# Plot
ggplot(sig_photo_dw_anno_summary, aes(x = as.factor(cluster), y = total_baseMean, fill = pathway)) +
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
# Join this dataframe with the original res_shrunk_dw_DF_anno to add the new column
res_shrunk_dw_DF_FULL <- res_shrunk_dw_DF_anno %>%
  rownames_to_column(var = "gene_id") %>%  # Convert rownames to a column for the join
  left_join(df_final_unique %>% select(Gene, IPRcode, Description, Model, ID, Function), by = c("gene_id" = "Gene")) %>%
  column_to_rownames(var = "gene_id") 

# Which of the genes are significantly DE? Now, I get the annotated version
sig_photo_dw_FULL <- res_shrunk_dw_DF_FULL %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(threshold == TRUE)  # This does the selecting for pvalue and log2fold change

write.csv(sig_photo_dw_FULL, file="df_sig_all_dw_ANNO_test3-2.csv")



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

# Join this dataframe with res_shrunk_dw_DF_FULL
res_shrunk_dw_DF_FULL_GO <- res_shrunk_dw_DF_FULL %>%
  rownames_to_column(var = "gene_id") %>%  # Convert rownames to a column for the join
  left_join(GO_df %>% select(gene_id, Term), by = c("gene_id" = "gene_id")) %>%
  column_to_rownames(var = "gene_id") 

View(res_shrunk_dw_DF_FULL_GO)

# Which of the genes are significantly DE? Now, I get the annotated version
sig_photo_dw_FULL_GO <- res_shrunk_dw_DF_FULL_GO %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(threshold == TRUE)  # This does the selecting for pvalue and log2fold change

write.csv(sig_photo_dw_FULL_GO, file="df_sig_all_dw_ANNO_test3-2.csv")


# Attach terms with most significant or most expressed list from before
# But I should really be doing this will all the sig expressed genes, and then match the clusters
#df_with_terms_subset <- df_with_terms[df_with_terms$gene_id %in% Top30_sig_list, ]
#write.csv(df_with_terms_subset, file="Top30_with_terms.csv")

# Assume deg_list is your list of differentially expressed genes
deg_list <- sig_res_dw$gene

# Prepare TERM2GENE and TERM2NAME
term2gene <- df_go_with_terms[, c("go_id", "gene_id")]
term2name <- bp_terms[, c("GOID", "TERM")]
colnames(term2name) <- c("ID", "Description")

# OR Prepare TERM2GENE and TERM2NAME for IPRcode
#term2gene <- df_ipr_unique[, c("go_id", "gene_id")]
#term2name <- df_ipr_unique[, c("go_id", "Description")]
#colnames(term2name) <- c("ID", "Description")


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

# Or do this

#all_genes_present <- all(deg_list %in% pep_ipr_only$Gene)

#if(all_genes_present) {
#  cat("All genes in deg_list are present in TERM2GENE.\n")
#} else {
#  cat("Some genes in deg_list are not present in TERM2GENE. Missing genes:\n")
#  print(setdiff(deg_list, pep_ipr_only$Gene))
#}


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

background_genes_df <- res_shrunk_dw %>%
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

write.csv(ego, file="df_dW_GO_BP_test3_2.csv")



## GSEA analyses
all_res_dw <- res_shrunk_dw %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

ranked_genes_df <- all_res_dw %>%
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

write.csv(gsea_result, file="df_dW_GSEA_BP_test3-2.csv")

# GSEA plot for a specific pathway
gseaplot(gsea_result, geneSetID = gsea_result$ID[1])




######## Analyzing expression levels (NOT Differential expression)
# create a vector of gene lengths 
myco_dw_raw <- mycounts_raw %>%
  dplyr::select(Length,
                S9A6PM,
                S9A6AM,
                S9A2AM,
                S9A10PM,
                S9A10AM,
                S6H6PM,
                S6H6AM,
                S6H2PM,
                S6H10PM,
                S6H10AM,
                S4J6PM,
                S4J6AM,
                S4J2PM,
                S4J2AM,
                S4J10PM,
                S4J10AM,
                S4H6PM,
                S4H2PM,
                S4H10PM,
                S4H10AM,
                S4B6PM,
                S4B6AM,
                S4B2PM,
                S4B2AM,
                S4B10PM,
                S4B10AM,
                S12N6PM,
                S12N6AM,
                S12N2PM,
                S12N2AM,
                S12N10PM,
                S12N10AM,
                S11L6PM,
                S11L6AM,
                S11L2PM,
                S11L2AM,
                S11L10PM,
                S11L10AM,
                S10F6PM,
                S10F6AM,
                S10F2PM,
                S10F2AM,
                S10F10AM)
View(myco_dw_raw)


######## Calculate TPM
# Get Gene length
geneLengths <- as.vector(subset(myco_dw_raw, select = c(Length)))

# Divide each gene count per sample by gene length, normalize by 1 million, normalize by total counts per sample 
x <- myco_dw_raw / geneLengths
myco_dw_tpm <- t( t(x) * 1e6 / colSums(x) )

View(myco_dw_tpm) # Does "Length have the same number?

## OR (same thing as above but easier to see what's going on)
#Calculate RPK
#gene_length <- as.numeric(as.vector(myco_dw_raw$Length))
#rpk <- myco_dw_raw / (gene_length / 1000)

#Calculate the scaling factor
#scaling_factor <- colSums(rpk) / 1e6

# Calculate TPM
#tpm <- sweep(rpk, 2, scaling_factor, FUN = "/")
#View(tpm)

# Delete the length column
myco_dw_tpm <- myco_dw_tpm[ ,-1]


#### Script to convert Time to Zeitgeber time
convert_to_zeitgeber <- function(time_str) {
  time_map <- c("10AM" = 2, "2PM" = 6, "6PM" = 10, "10PM" = 14, "2AM" = 18, "6AM" = 22)
  return(time_map[time_str])
}


#### Subset PPC genes
ppc_genes <- c("g7768", "g11025", "g9957", "g22069")

# Load cleaned version
ortho_long <- read.csv(file = 'photosynthetic_genes_orthologs_longer.csv')

ortho_long_subset <- ortho_long[ortho_long$Cislong_ortho %in% ppc_genes, ]
View(ortho_long_subset)

# Change WD
setwd("~/Dropbox/Yale/Research/RNAseq/TPM_graphs_DW_NEW")

# Step 1: Extract count data for the specific gene
gene_id <- "g11025"  # Replace with your gene of interest
gene_df <- myco_dw_tpm[rownames(myco_dw_tpm) == gene_id, , drop = FALSE]

# Ensure the gene exists in the dataset
if (nrow(gene_df) == 0) {
  stop(paste("Gene ID", gene_id, "not found in myco_dw_tpm."))
}

# Step 2: Reshape the data and change the first column name
gene_df <- as.data.frame(t(gene_df))
gene_df$SampleID <- rownames(gene_df)

colnames(gene_df)[1] = "counts"

# Step 3: Merge the count data with the metadata-
# Now,`merged_data` has the `count`, `Time`, and `Treatment` columns along with other metadata
merged_data <- merge(gene_df, meta_dw, by = "SampleID")

# Step 4: Convert Time to Zeitgeber_Time
merged_data$Zeitgeber_Time <- sapply(merged_data$Time, convert_to_zeitgeber)

# Step 5: Make Summary of the data to compute means and SD
summary_df <- merged_data %>%
  group_by(Zeitgeber_Time, Treatment) %>%
  summarise(mean_counts = mean(counts), se_counts = sd(counts)/sqrt(n()))

# Step 6: Plotting
p <- ggplot(summary_df, aes(x = Zeitgeber_Time, y = mean_counts, color = Treatment, group = Treatment)) + 
  geom_errorbar(aes(ymin = mean_counts - se_counts, ymax = mean_counts + se_counts), width = 0.5, linewidth = 2) +
  geom_line(linewidth = 2) +  
  scale_color_manual(labels = c("Watered", "Drought"), values=c("#619CFF", "#F8766D")) +
  ylab("TPM") +
  xlab("Zeitgeber Time (hours)") +
  scale_x_continuous(breaks = seq(4, 24, by = 4), limits = c(1, 24)) +
  annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf, fill = "lightgrey", alpha = 0.3) +
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

#### Loop through each gene in the subset or the full dataset

# Change WD
setwd("~/Dropbox/Yale/Research/RNAseq/TPM_graphs_DW_NEW")

# Loop
for (i in 1:nrow(ortho_long)) {
  gene_id <- ortho_long$Cislong_ortho[i]
  label <- ortho_long$Label[i]
  
  # Step 1: Extract count data for the specific gene
  gene_df <- myco_dw_tpm[rownames(myco_dw_tpm) == gene_id, , drop = FALSE]
  
  # Ensure the gene exists in the dataset
  if (nrow(gene_df) == 0) {
    warning(paste("Gene ID", gene_id, "not found in myco_dw_tpm. Skipping..."))
    next
  }
  
  # Step 2: Reshape the data using pivot_longer
  gene_df <- as.data.frame(t(gene_df))
  gene_df$SampleID <- rownames(gene_df)
  
  colnames(gene_df)[1] = "counts"
  
  # Step 3: Merge the count data with the metadata
  merged_data <- merge(gene_df, meta_dw, by = "SampleID")
  
  # Convert Time to Zeitgeber time
  merged_data$Zeitgeber_Time <- sapply(merged_data$Time, convert_to_zeitgeber)
  
  # Summarize data by Zeitgeber time and treatment
  summary_df <- merged_data %>%
    group_by(Zeitgeber_Time, Treatment) %>%
    summarise(mean_counts = mean(counts), se_counts = sd(counts)/sqrt(n()))
  
  # Step 4: Plotting
  p <- ggplot(summary_df, aes(x = Zeitgeber_Time, y = mean_counts, color = Treatment, group = Treatment)) + 
    geom_errorbar(aes(ymin = mean_counts - se_counts, ymax = mean_counts + se_counts), width = 0.5, linewidth = 2.5) +
    geom_line(linewidth = 2.5) +  
    scale_color_manual(labels = c("Watered", "Drought"), values=c("#619CFF", "#F8766D")) +
    ylab("TPM") +
    xlab("Zeitgeber time (hours)") +
    scale_x_continuous(breaks = seq(4, 24, by = 4), limits = c(1, 24)) +
    annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf, fill = "lightgrey", alpha = 0.3) +
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




## fin