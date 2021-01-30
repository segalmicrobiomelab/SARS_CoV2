library(ggplot2)
library(readxl)
library(edgeR)
library(Seurat)
library(xCell)
library(reshape2)
library(ggpubr)
library(scales)

setwd("~/Documents/SARS-Cov2_metatranscriptome_project/Deconvolution/")
data <- read.delim("RNA_Host_Transcriptome_ALL.txt", stringsAsFactors=FALSE, row.names=1)
meta <- read_xlsx("../files_from_leo/201216.RNA.xlsx")

#### Prepare single-cell reference signature matrix for CIBERSORT ####
## Use study Liao et al (2020) Nature Medicine ##
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926
# https://raw.githubusercontent.com/zhangzlab/covid_balf/master/all.cell.annotation.meta.txt

# A more detailed annotation #
ref_h5_files <- list.files("Liao2020NatureMedicine_single-cell-ref/GSE145926_RAW/", pattern=".h5")
scref <- list()
for (i in 1:length(ref_h5_files)) {
  scref[[i]] <- Read10X_h5(paste("Liao2020NatureMedicine_single-cell-ref/GSE145926_RAW/", ref_h5_files[i], sep=""))
  names(scref)[i] <- gsub("_filtered_feature_bc_matrix.h5", "", ref_h5_files[i])
}

scref_celltype <- read.delim("Liao2020NatureMedicine_single-cell-ref/all.cell.annotation.meta.txt", stringsAsFactors=FALSE)
scref_celltype_myeloid <- read.delim("Liao2020NatureMedicine_single-cell-ref/myeloid.cell.annotation.meta.txt", stringsAsFactors=FALSE)
scref_celltype_NKT <- read.delim("Liao2020NatureMedicine_single-cell-ref/NKT.cell.annotation.meta.txt", stringsAsFactors=FALSE)
# reformulate cell-type meta data
scref_celltype_other <- scref_celltype[!scref_celltype$ID %in% c(scref_celltype_myeloid$ID, scref_celltype_NKT$ID), ]
scref_celltype_old <- scref_celltype
scref_celltype <- rbind(scref_celltype_myeloid, scref_celltype_NKT, scref_celltype_other)
unique(scref_celltype$celltype)
# Change macrophage groups into actual names
scref_celltype$celltype <- gsub("Group1|Group2", "Macrophages M1", scref_celltype$celltype)
scref_celltype$celltype <- gsub("Group3", "Macrophages M2", scref_celltype$celltype)
scref_celltype$celltype <- gsub("Group4", "Alveolar macrophages", scref_celltype$celltype)
# Remove "Doublets" and "Uncertain" T-cells
scref_celltype <- scref_celltype[!scref_celltype$celltype %in% c("Uncertain", "Doublets"), ]

# For each sample and each cell type -> sample 10 cells
scref_subset <- list()
for (i in 1:length(scref)) {
  i_sample_fullID <- names(scref)[[i]]
  i_sample <- strsplit(names(scref)[[i]], split="_")[[1]][2]
  i_sample_data <- as.data.frame(scref[[i_sample_fullID]])
  
  i_df <- scref_celltype[scref_celltype$sample == i_sample, ]
  i_output <- data.frame()
  i_celltypes <- character()
  for (j in 1:length(unique(as.character(i_df$celltype)))) {
    ij_celltype <- unique(as.character(i_df$celltype))[j]
    if (nrow(i_df[i_df$celltype == ij_celltype, ]) < 10) { # If cell number < 10, then take all
      size <- nrow(i_df[i_df$celltype == ij_celltype, ])
    } else {
      size <- 10
    }
    ij_cellID <- sample(i_df[i_df$celltype == ij_celltype, "ID"], size=size, replace=FALSE)
    ij_cellID <- paste(unlist(strsplit(ij_cellID, split="_"))[seq(1, length(unlist(strsplit(ij_cellID, split="_"))), 2)], "-1", sep="")
    ij_sample_data <- as.data.frame(i_sample_data[, colnames(i_sample_data) %in% ij_cellID])
    if (ncol(ij_sample_data) == 1) {
      rownames(ij_sample_data) <- rownames(i_sample_data)
      colnames(ij_sample_data) <- ij_cellID
    }
    ij_sample_data$Gene <- rownames(ij_sample_data)
    if (ncol(i_output) == 0) {
      i_output <- ij_sample_data
    } else {
      i_output <- merge(i_output, ij_sample_data, by="Gene", all=TRUE)
    }
    i_celltypes <- c(i_celltypes, rep(ij_celltype, size))
  }
  i_output <- rbind(i_output, c("CellType", rep(i_celltypes, ncol(i_output)-1)))
  scref_subset[[i_sample]] <- i_output
} 

scref_subset_df <- data.frame()
for (i in 1:length(scref_subset)) {
  i_df <- scref_subset[[i]]
  i_df <- i_df[i_df$Gene != "nCoV", ] # remove nCoV counts
  if (i == 1) {
    scref_subset_df <- i_df
  } else {
    scref_subset_df <- merge(scref_subset_df, i_df, by="Gene", all=TRUE)
  }
}
ct_idx <- which(scref_subset_df$Gene == "CellType")
scref_subset_df <- scref_subset_df[c(ct_idx, 1:(ct_idx-1), (ct_idx+1):nrow(scref_subset_df)), ]
scref_subset_df[1,1] <- "GeneSymbol"
write.table(scref_subset_df, "Liao2020NatureMedicine_single-cell-ref/Liao2020_detailed-annot_scref_matrix.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
table(as.character(scref_subset_df[1, ])) # T cells have the minimum cell number as 20



#### Prepare gene lenght information ####
# Export entrex gene ID for convertion #
load("hg19_geneID_convertion/COVIDRNAseqCounts.Rdata")
genelist <- mycounts2[["annotation"]]
write.table(as.data.frame(genelist[, "GeneID"]), file="hg19_entrezID.txt", 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# Reload converted entrex gene ID to gene symbols #
DAVID_genesymbol <- read.delim("hg19_geneID_convertion/DAVID_convertion.txt", header=TRUE, stringsAsFactors=FALSE)
data <- read.delim("RNA_Host_Transcriptome_ALL.txt", stringsAsFactors=FALSE)

sum(!is.na(DAVID_genesymbol$To)) 
sum(! DAVID_genesymbol$To %in% data$X) 

# Export genesymbol and length information #
gene_symbol_length <- merge(genelist[, c("GeneID", "Length")], DAVID_genesymbol[, c("From", "To")], by.x="GeneID", by.y="From")
colnames(gene_symbol_length) <- c("Geneid", "Length", "Genename")
gene_symbol_length <- gene_symbol_length[, c(1,3,2)]
gene_symbol_length$Genename <- toupper(gene_symbol_length$Genename)
length(unique(gene_symbol_length$Genename)) 
write.table(gene_symbol_length, file="gene_symbol_length_DAVID-convertion.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

rm(DAVID_genesymbol, mycounts2)



#### Prepare mixture input data TPM - BAL samples only and after moving of a few samples ####
rownames(data) <- data$X
data$X <- NULL
## Only Keep the BAL samples ##
data <- data[, grepl("BAL", colnames(data))]
## Remove a few samples according to the metadata ##
sample_keep <- as.character(as.data.frame(meta[meta$Final_Paper == 1, "Study_Linked.ID"])$Study_Linked.ID)
data <- data[, colnames(data) %in% sample_keep]
sum(as.data.frame(meta[meta$Sample.Type == "BAL" & meta$Final_Paper == 1, "Study_Linked.ID"])$Study_Linked.ID %in% colnames(data)) == ncol(data) # TRUE

## Filter and normalize gene counts ##
# CPM filtering
plot(log10(colSums(data)))
mean(colSums(data)) 
median(colSums(data)) 
data_cpm <- cpm(data) 
data_keep <- rowSums(data_cpm >= 1) >= 2 # at least 1cpm in at least 2 samples
data <- data[data_keep, ]
dim(data) 

# Normalization by gene length (using TPM)
data_lennorm <- data
data_lennorm$name <- row.names(data_lennorm)
data_lennorm <- merge(gene_symbol_length[, c("Genename", "Length")], data_lennorm, 
                      by.x="Genename", by.y="name")

for (i in 3:ncol(data_lennorm)) {
  data_lennorm[, i] <- data_lennorm[, i] / (data_lennorm$Length/1000)
}

length(unique(data_lennorm$Genename)) == nrow(data_lennorm) # TRUE
rownames(data_lennorm) <- data_lennorm$Genename
data_lennorm <- data_lennorm[, 3:ncol(data_lennorm)]
data_tpm <- cpm(data_lennorm)
colSums(data_tpm) #1e6
dim(data_tpm) # 17340 genes for 118 samples 

data_tpm <- as.data.frame(data_tpm)
data_tpm$Gene <- rownames(data_tpm)
data_tpm <- data_tpm[, c(ncol(data_tpm), 1:(ncol(data_tpm)-1))]
write.table(data_tpm, file="RNA_Host_Transcriptome_BAL_TPM_Final.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)



#### xCell ####
data_tpm <- read.delim("RNA_Host_Transcriptome_BAL_TPM_Final.txt", stringsAsFactors=FALSE)
rownames(data_tpm) <- data_tpm$Gene
data_tpm$Gene <- NULL

cell.types.use <- c("Epithelial cells", 
                    "Macrophages", "Macrophages M1", "Macrophages M2",
                    "Neutrophils", "Mast cells", 
                    "DC", "cDC", "pDC", "NK cells", 
                    "B-cells", "Memory B-cells", "Plasma cells", 
                    "CD4+ T-cells", "Th1 cells", "Th2 cells", "Tregs", "Tgd cells",
                    "CD8+ T-cells", "NKT")
xcell_scores <- xCellAnalysis(as.matrix(data_tpm), cell.types.use=cell.types.use, rnaseq=TRUE)

meta <- read_xlsx("../files_from_leo/201216.RNA.xlsx")
meta <- meta[, c("Study_Linked.ID", "Sample.Type", "Final_Paper", "@3_groups")] 
colnames(meta) <- c("Study_Linked.ID", "Sample.Type", "Final_Paper", "groups")
meta$groups <- gsub("≤.28.day.on.vent", "<28.day.on.vent", meta$groups)

# dot and box plot
xcell_scores_melt <- melt(xcell_scores)
colnames(xcell_scores_melt) <- c("celltype", "Sample.ID", "score")
xcell_scores_melt <- merge(meta, xcell_scores_melt, by.x="Study_Linked.ID", by.y="Sample.ID", all.y=TRUE)
xcell_scores_melt$celltype <- factor(xcell_scores_melt$celltype, levels=cell.types.use)
xcell_scores_melt$groups <- factor(xcell_scores_melt$groups, levels=c("<28.day.on.vent", ">28.day.on.vent", "Dead"))

ggplot(xcell_scores_melt, aes(x=groups, y=score)) +
  geom_jitter(aes(color=groups), size=2) +
  geom_boxplot(fill="transparent", outlier.colour="transparent") +
  facet_wrap(~celltype, ncol=5, scale="free_y") +
  scale_color_manual(values=c("#41a01b", "#fecb66", "#c4057b")) +
  theme(panel.background=element_rect(fill="transparent", color="grey90"),
        panel.grid.major=element_line(color="grey90"),
        text=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="Group", y="xCell enrichment socre", color="Group") +
  scale_y_continuous(labels=scientific)

# Statistical tests - wilcoxon and BH correction
xcell_wilcx <- data.frame(matrix(ncol=3, nrow=0), stringsAsFactors=FALSE)
for (i in 1:length(cell.types.use)) {
  xcell_score_i <- xcell_scores_melt[as.character(xcell_scores_melt$celltype) == cell.types.use[i], ]
  pairwise_p_i <- melt(pairwise.wilcox.test(xcell_score_i$score, xcell_score_i$groups, p.adjust.method="none")$p.value)
  pairwise_p_i <- pairwise_p_i[! is.na(pairwise_p_i$value), ]
  pairwise_p_i$comparison <- paste(pairwise_p_i$Var1, pairwise_p_i$Var2, sep="VS")
  pairwise_p_i$celltype <- cell.types.use[i]
  pairwise_p_i <- pairwise_p_i[, c(5,4,3)]
  colnames(pairwise_p_i)[3] <- "pval"
  xcell_wilcx <- rbind(xcell_wilcx, pairwise_p_i)
}
xcell_wilcx$pval_adj <- p.adjust(xcell_wilcx$pval, method="BH")
xcell_wilcx_sig <- xcell_wilcx[xcell_wilcx$pval_adj <= 0.05, ]

# Export 
write.table(xcell_scores, "118BALsamples_Final/xCell/xcell_scores_BALsamples.txt", quote=FALSE, row.names=TRUE, sep="\t")
write.table(xcell_wilcx, "118BALsamples_Final/xCell/xcell_score_wilcx_BALsamples.txt", quote=FALSE, row.names=FALSE, sep="\t")



#### CIBERSORTx absolute mode results ####
cibersort_absolute_res <- read.delim("118BALsamples_Final/CIBERSORTx/CIBERSORTx_Job10_Adjusted.txt", stringsAsFactors=FALSE)
# Check the p-value for deconvolution
plot(cibersort_absolute_res$P.value) # all 0

cibersort_absolute_res_new <- cibersort_absolute_res[, 1:19]
cibersort_absolute_res_new <- merge(meta, cibersort_absolute_res_new, by.x="Study_Linked.ID", by.y="Mixture")
cibersort_absolute_res_new <- melt(cibersort_absolute_res_new, id.vars=colnames(meta))
colnames(cibersort_absolute_res_new)[5:6] <- c("celltype", "frac")
cibersort_absolute_res_new$groups <- factor(cibersort_absolute_res_new$groups, levels=c("<28.day.on.vent", ">28.day.on.vent", "Dead"))
cibersort_absolute_res_new$celltype <- factor(cibersort_absolute_res_new$celltype, 
                                              levels=c("Epithelial", "Macrophages", "Macrophages.M1", "Macrophages.M2", "Alveolar.macrophages", 
                                                       "Neutrophil", "Mast", "mDC", "pDC", "NK", 
                                                       "B", "Plasma", "T", "innate.T", "CCR7..T", "Treg", "CD8.T", "Proliferating.T"))

ggplot(cibersort_absolute_res_new, aes(x=groups, y=frac)) +
  geom_jitter(aes(color=groups), size=2) +
  geom_boxplot(fill="transparent", outlier.colour="transparent") +
  facet_wrap(~celltype, ncol=5, scale="free_y") +
  scale_color_manual(values=c("#41a01b", "#fecb66", "#c4057b")) +
  theme(panel.background=element_rect(fill="transparent", color="grey90"),
        panel.grid.major=element_line(color="grey90"),
        text=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="Group", y="CIBERSORTx absolute cell fraction score", color="Group")

# Statistical tests - wilcoxon and BH correction
cibersort_cell.types.use <- levels(cibersort_absolute_res_new$celltype)
# test for normality
shapiro.test(cibersort_absolute_res_new$frac) #p-value < 2.2e-16
ggqqplot(cibersort_absolute_res_new$frac)

cibersort_absolute_wilcx <- data.frame(matrix(ncol=4, nrow=0), stringsAsFactors=FALSE)
for (i in 1:length(cibersort_cell.types.use)) {
  cibersortx_frac_i <- cibersort_absolute_res_new[as.character(cibersort_absolute_res_new$celltype) == cibersort_cell.types.use[i], ]
  pairwise_p_i <- melt(pairwise.wilcox.test(cibersortx_frac_i$frac, cibersortx_frac_i$groups, p.adjust.method="none")$p.value)
  pairwise_p_i <- pairwise_p_i[! is.na(pairwise_p_i$value), ]
  pairwise_p_i$comparison <- paste(pairwise_p_i$Var1, pairwise_p_i$Var2, sep="VS")
  pairwise_p_i$celltype <- cibersort_cell.types.use[i]
  pairwise_p_i <- pairwise_p_i[, c(5,4,3)]
  colnames(pairwise_p_i)[3] <- "pval"
  cibersort_absolute_wilcx <- rbind(cibersort_absolute_wilcx, pairwise_p_i)
}
cibersort_absolute_wilcx$pval_adj <- p.adjust(cibersort_absolute_wilcx$pval, method="BH")
cibersort_absolute_wilcx_sig <- cibersort_absolute_wilcx[cibersort_absolute_wilcx$pval_adj <= 0.05, ]

# Export 
write.table(cibersort_absolute_wilcx, "118BALsamples_Final/CIBERSORTx/cibersort_LiaoRef_Job10_absolute_wilcx_BALsamples.txt", quote=FALSE, row.names=FALSE, sep="\t")



#### Kruskal-Walis rank sum test for each cell type between all the groups ####
## Load data ##
meta <- read_xlsx("../files_from_leo/201216.RNA.xlsx")
meta <- meta[, c("Study_Linked.ID", "Sample.Type", "Final_Paper", "@3_groups")] 
colnames(meta) <- c("Study_Linked.ID", "Sample.Type", "Final_Paper", "groups")
meta$groups <- gsub("≤.28.day.on.vent", "<28.day.on.vent", meta$groups)

xcell_scores <- read.delim("118BALsamples_Final/xCell/xcell_scores_BALsamples.txt", stringsAsFactors=FALSE)
xcell_scores <- as.data.frame(t(xcell_scores))
xcell_scores$Mixture <- rownames(xcell_scores)
cibersort_absolute_res <- read.delim("118BALsamples_Final/CIBERSORTx/CIBERSORTx_Job10_Adjusted.txt", stringsAsFactors=FALSE)
cibersort_absolute_res <- cibersort_absolute_res[, 1:19]
colnames(cibersort_absolute_res) <- c("Mixture",
                                      "Macrophages M1", "Alveolar macrophages", "Macrophages M2",
                                      "Innate T-cells", "Tregs", "CD8+ T-cells", "NK cells",
                                      "CCR7+ T-cells", "Proliferating T-cells", "Other macrophages",
                                      "Epithelial cells", "B-cells", "cDC", "pDC", "Neutrophils",
                                      "Plasma cells", "Mast cells", "Other T-cells")

## Kruskal-Walis test ##
kruskal_test_fun <- function(input) {
  input_melt <- melt(input, id.vars="Mixture")
  colnames(input_melt) <- c("Study_Linked.ID", "celltype", "score")
  input_melt <- merge(meta, input_melt, by="Study_Linked.ID")
  input_melt$groups <- factor(input_melt$groups, levels=c("<28.day.on.vent", ">28.day.on.vent", "Dead"))
  
  output <- data.frame(matrix(ncol=2, nrow=0), stringsAsFactors=FALSE)
  for (i in 1:length(unique(as.character(input_melt$celltype)))) {
    score_i <- input_melt[as.character(input_melt$celltype) == unique(as.character(input_melt$celltype))[i], ]
    pairwise_p_i <- kruskal.test(score_i$score, g=score_i$groups)$p.value
    output[i, ] <- c(unique(as.character(score_i$celltype)), pairwise_p_i)
  }
  colnames(output) <- c("celltype", "pval")
  output$adj_pval <- p.adjust(output$pval, method="BH")
  return(output)
}

xcell_kruskal <- kruskal_test_fun(xcell_scores)
cibersortxabs_kruskal <- kruskal_test_fun(cibersort_absolute_res)

cibersortxabs_xcell_kruskal <- merge(cibersortxabs_kruskal, xcell_kruskal, 
                                     by="celltype", suffixes=c(".CIBERSORTxabs", ".xCell"), all=TRUE)

## Export ##
write.table(xcell_kruskal, "118BALsamples_Final/xCell/xCell_outcome-group_kruskal-Walis-test.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(cibersortxabs_kruskal, "118BALsamples_Final/CIBERSORTx/CIBERSORTxabs_outcome-group_kruskal-Walis-test.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(cibersortxabs_xcell_kruskal, "118BALsamples_Final/CIBERSORTxabs-xCell_outcome-group_kruskal-Walis-test_summary.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)



#### re-plot CIBERSORTx absolute mode results for Macrophages M1, Mast cells, innate T cells and CCR7 ####
meta <- read_xlsx("../files_from_leo/201216.RNA.xlsx")
meta <- meta[, c("Study_Linked.ID", "Sample.Type", "Final_Paper", "@3_groups")] 
colnames(meta) <- c("Study_Linked.ID", "Sample.Type", "Final_Paper", "groups")
meta$groups <- gsub("≤.28.day.on.vent", "<28.day.on.vent", meta$groups)

cibersort_absolute_res <- read.delim("118BALsamples_Final/CIBERSORTx/CIBERSORTx_Job10_Adjusted.txt", stringsAsFactors=FALSE)
# Check the p-value for deconvolution
plot(cibersort_absolute_res$P.value) # all 0

cibersort_absolute_res_new <- cibersort_absolute_res[, c("Mixture", "Macrophages.M1", "Mast", "innate.T", "CCR7..T")]
colnames(cibersort_absolute_res_new) <- c("Mixture", "Macrophages M1", "Mast cells", "Innate T-cells", "CCR7+ T-cells")
cibersort_absolute_res_new <- merge(meta, cibersort_absolute_res_new, by.x="Study_Linked.ID", by.y="Mixture")
cibersort_absolute_res_new <- melt(cibersort_absolute_res_new, id.vars=colnames(meta))
colnames(cibersort_absolute_res_new)[5:6] <- c("celltype", "frac")
cibersort_absolute_res_new$groups <- factor(cibersort_absolute_res_new$groups, levels=c("<28.day.on.vent", ">28.day.on.vent", "Dead"))
cibersort_absolute_res_new$celltype <- factor(cibersort_absolute_res_new$celltype, 
                                              levels=c("Macrophages M1", "Mast cells", "Innate T-cells", "CCR7+ T-cells"))
                                              
ggplot(cibersort_absolute_res_new, aes(x=groups, y=frac)) +
  geom_jitter(aes(color=groups), size=3) +
  geom_boxplot(fill="transparent", outlier.colour="transparent") +
  facet_wrap(~celltype, ncol=4, scale="free_y") +
  scale_color_manual(values=c("#41a01b", "#fecb66", "#c4057b")) +
  theme(panel.background=element_rect(fill="transparent", color="grey90"),
        panel.grid.major=element_line(color="grey90"),
        text=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position="none") +
  labs(x="Group", y="CIBERSORTx absolute cell fraction score", color="Group")
# size 10x5
