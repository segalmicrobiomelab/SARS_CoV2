#Load Packages
library(DESeq2)
library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(pathfindR)
library(scales)
library(data.table)
library(fBasics)
library(forcats)
library(omu)
library(maptools)
library(phyloseq)
library(vegan)
#library(ggpmisc)
library(dplyr)
library(tibble)
library(formattable)
library("htmltools")
library("webshot")    
library(splitstackshape)
library(decontam)
library(dplyr)
library(grid)
library(cowplot)
library(wesanderson)
library(colorspace)
#library(egg)


#Set Theme for Figures
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")

#Choose Alpha/FDR
alpha = 0.05

#Create Function For Analysis
deseq <- function(metadata,counts,sampletype,comparison) {
        #Load MetaData
        coldata <- read.delim2(metadata, sep="\t")
        rownames(coldata) <- coldata$Study_Linked_ID
        #Fix Three Groups
        coldata$three_groups <- ifelse(coldata$X.3_groups=="\xb2.28.day.on.vent","Less_Than_28_days_on_vent",
                                ifelse(coldata$X.3_groups==">28.day.on.vent","Greater_Than_28_days_on_vent",
                                as.character(coldata$X.3_groups)))
        #load Count Data
        mycounts <-read.delim2(paste0(counts,".txt"), sep="\t", row.names=1)
        #Find matching sample ID for both datasets
        needed<-which(rownames(coldata) %in% colnames(mycounts))    
        #keep only matching IDs from count data
        coldata2<-coldata[needed,]
        #Order Meta Data by SampleId
        coldata2 <- coldata2[order(coldata2$Study_Linked_ID),]
        #keep only matching IDs from count data
        wanted<-which(colnames(mycounts) %in% rownames(coldata))    
        #keep only matching IDs from count data
        mycounts2<-mycounts[,wanted]
        #Order Count Data by SampleID
        mycounts2 <-mycounts2[, order(colnames(mycounts2))]
        #Convert any NAs to 0
        mycounts2[is.na(mycounts2)] <- 0
        #Copy of Count Table
        mycounts3 <- mycounts2
        #Convert Count Table into a Numeic Data Frame
        d1 = data.frame(lapply(mycounts3, function(x) as.numeric(as.character(x))),
                           check.names=F, row.names = rownames(mycounts3))
        #Convert Data to Integers to Run DESEq
        d1[] <- lapply(d1, as.integer)
        ddsv <- DESeqDataSetFromMatrix(countData = d1,
                              colData = coldata2,
                              design= ~ Sample.Type)
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsv), 1, gm_mean)
        #Estimate Factors of DESeq Object
        ddsv <- estimateSizeFactors(ddsv, geoMeans = geoMeans)
        #Transforming data - Option #2 is variance stabilizing transformation
        vsdv <- varianceStabilizingTransformation(ddsv)
        #----------------------
        ##Create graph of Reads
        #----------------------
        #Sum number or reads per sample
        summary <- as.data.frame(rowSums(t(assay(ddsv))))
        #Merge Reads Data with MetaData
        require(data.table)
        Reads <- as.data.frame(merge(x = summary, y = colData(ddsv), by = "row.names", all.x = TRUE))
        #Rename Column of Reads
        colnames(Reads)[colnames(Reads)=="rowSums.t.assay.ddsv..."] <- "Reads"
        #Set Order Of Figure
        Reads$or <-ifelse(Reads$Sample.Type=="BKG", 1,NA)
        Reads$or <-ifelse(Reads$Sample.Type=="BAL",2 ,Reads$or)
        Reads$or <-ifelse(Reads$Sample.Type=="UA",3 ,Reads$or)
        #Create Figure
            ggsave(filename=paste0(counts,".Reads.pdf"),
            ggplot(Reads, aes(x= reorder(Sample.Type, +or), y=Reads, fill=Sample.Type)) + 
            stat_boxplot(geom ='errorbar', width=0.1)+
            geom_boxplot(outlier.shape = NA, width=0.5)+
            geom_jitter(shape=1, position=position_jitter(0.2))+
            scale_fill_manual(values=c("#1700F5", "#BEBEBE","#F7A501")) +
            scale_x_discrete(labels = c('BKG','BAL','UA'))+ 
            scale_y_continuous(name="Reads",trans="log10", breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
            xlab("Sample Type")+
            theme,
            height = 7, width = 5)
        #----------------------
        ##PCOA PLOT
        #----------------------
        #Create Distance Matrix
        vsdv0 <- ifelse(assay(vsdv)<0,0,assay(vsdv))
        vegdist   = vegdist(t(vsdv0), method="bray")
        #Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
        CmdScale <- cmdscale(vegdist, k =10)
        #calculated Sample variance for each PC
        vars <- apply(CmdScale, 2, var)
        #Create Variable with the Percent Variance
        percentVar <- round(100 * (vars/sum(vars)))
        #Merge PC Data with MetaData
        require(data.table)
        newResults <- merge(x = CmdScale, y = colData(vsdv), by = "row.names", all.x = TRUE)
        #Rename Variables for PC1 and PC2
        colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
        colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
        colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
        #Calculate the Centroid Value
        centroids <- aggregate(cbind(PC1,PC2)~ Sample.Type,data= newResults, mean)
        #Merge the Centroid Data into the PCOA Data
        newResults <- merge(newResults,centroids,by="Sample.Type",suffixes=c("",".centroid"))
        #Create Data For Statistics
        data.adonis <- data.frame(colData(vsdv))
        #Run the Statistics
        samplepermanova <- adonis(vegdist ~ Sample.Type, data.adonis)
        samplepermanova <- as.data.frame(samplepermanova$aov.tab)
        samplepermanova <- samplepermanova$'Pr(>F)'[1]
        #PLOT IT
            ggsave(filename=paste0(counts,".BRAY_vsd_sampletype_PERMANOVA_",samplepermanova,".pdf"),
            ggplot(newResults, aes(PC1, PC2,color=Sample.Type)) + # Graph PC1 and PC2
            xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
            ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
            scale_color_manual(values=c("#1700F5", "#BEBEBE", "#F7A501")) + 
            geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Sample.Type))+ 
            geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=Sample.Type), size=10) + #BKG DFW Lung.Tissue.In Lung.Tissue.UnIn Mock
            geom_point(data=newResults,aes(color=Sample.Type),size=5,alpha=0.5) + # Set the size of the points
            theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
            panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
            panel.grid.minor = element_blank(),strip.background=element_blank(),
            axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
            axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),
            plot.margin=unit(c(1,1,1,1),"line"),legend.position="none"),
            height = 10, width = 10)
        #----------------------
        ##Alpha Diversity
        #----------------------
        #Calcultes Shannon Diversity
        ddsv$Shannon = diversity(vsdv0, index = "shannon", MARGIN = 2, base = exp(1))
        #Convert to data frame for ggplot
        shannon = as.data.frame(colData(ddsv))
        #Remove any zero values
        shannon[shannon==0] <- NA
        #Set Order Of Figure
        shannon$or <-ifelse(shannon$Sample.Type=="BKG", 1,NA)
        shannon$or <-ifelse(shannon$Sample.Type=="BAL",2 ,shannon$or)
        shannon$or <-ifelse(shannon$Sample.Type=="UA",3 ,shannon$or)
        #Check Statistics
        sampleshannon <- kruskal.test(Shannon ~ Sample.Type, data = shannon)
        sampleshannon <- sampleshannon$p.value
        #Plot It
            ggsave(filename=paste0(counts,".SHANNON_vsd_sampletype_Pvalue_",sampleshannon,".pdf"),
            ggplot(shannon, aes(x= reorder(Sample.Type, +or), y=Shannon, fill=Sample.Type)) + 
            stat_boxplot(geom ='errorbar', width=0.1)+
            geom_boxplot(outlier.shape = NA, width=0.5)+
            #geom_boxplot(aes(ymin=..lower.., ymax=..upper..))+
            geom_jitter(shape=1, position=position_jitter(0.2))+
            scale_fill_manual(values=c("#1700F5", "#BEBEBE", "#F7A501")) +
            scale_x_discrete(labels = c('BKG','BAL','UA'))+ 
            #geom_text_repel(aes(label=ifelse(res$sig > 3 , as.character(res$o),'')),size=3,force=25) +
            #stat_summary(fun.y = max, colour = "red", geom = "point", size = 2)+
            #stat_summary(aes(label=otu), fun.y=max, geom="text", size=6, hjust = -0.3)+
            ylab("Shannon Diversity") + 
            xlab("Sample Type")+
            #geom_text(data=subset(res, sig==max(sig)), aes(label=otu)) +
            #geom_point(color=cols) +
            #geom_text(aes(x=sample, y=sig, label = o), data2,  col = 'red') +
            theme, 
            height = 7, width = 5)
}

