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
        #----------------------
        ##DECONTAM
        #----------------------
        GenusData <-as.data.frame(assay(ddsv)) #pruned to selected Genuses based on abundance        
        #Create Relative Abundance Table
        df <-
            GenusData %>% 
                rownames_to_column('gs') %>%
                group_by(gs) %>% 
                summarise_all(funs(sum)) %>%
                mutate_if(is.numeric, funs(./sum(.))) %>%
                column_to_rownames('gs')
        #Create a TRUE/FALSE Variable for Subject and Control
        coldata2$is.neg <- coldata2$Sample.Type=="BKG"
        #Transpose to matrix
        otus <- t(as.matrix(df))
        ## extract only those samples in common between the two tables
        common.sample.ids <- intersect(rownames(coldata2), rownames(otus))
        otus <- otus[common.sample.ids,]
        coldata2 <- coldata2[common.sample.ids,]
        #Check Contaminants based on BKG samples
        contamdf.prev <- isContaminant(otus, method="prevalence", neg=coldata2$is.neg,threshold=0.5)
        table(contamdf.prev$contaminant)
        #Select only BKG Samples
        BKG <- df %>% select(rownames(coldata2[coldata2$Sample.Type=='BKG',]))
        #List of the Contaminants Identified
        otu.to.save <-as.character(rownames(contamdf.prev[contamdf.prev$contaminant==TRUE,]))
        #sum of BKG Abundances
        BKG <-rowSums(BKG)
        #Only select contaminants
        BKG <- BKG[otu.to.save]
        #Sort table by Total Abundance
        BKG <- sort(BKG,decreasing=TRUE)
        #transform to Data Frame
        BKG <- as.data.frame(BKG)
        #Select top 10
        BKG <- BKG %>% dplyr::slice(1:10)
        #write Table
        write.table(BKG,file=paste0(counts,".top10_contaminants.txt"), sep="\t", col.names = NA, row.names = TRUE)
        #----------------------
        ##DESEQ
        #----------------------
        #Create Deseq object for BAL analysis
        ddsvbal <- DESeqDataSetFromMatrix(countData = d1,
                              colData = coldata2,
                              design= ~ Sample.Type)
        #Subset just the BAL
        ddsvbal <- ddsvbal[, ddsvbal$Sample.Type %in% "BAL"]
        #Covert Variable to Factor
        ddsvbal$three_groups <- as.factor(ddsvbal$three_groups)
        #Create Seperate tables with each of the different compairsons
        ddsvbal1 <- ddsvbal[, ddsvbal$three_groups %in% c("Dead","Less_Than_28_days_on_vent")]
        ddsvbal2 <- ddsvbal[, ddsvbal$three_groups %in% c("Less_Than_28_days_on_vent","Greater_Than_28_days_on_vent")]
        ddsvbal3 <- ddsvbal[, ddsvbal$three_groups %in% c("Greater_Than_28_days_on_vent","Dead")]
        #Create Tables of First Comparison
        balcounts1 <- assay(ddsvbal1)
        balmeta1   <- colData(ddsvbal1)
        #Create Tables of Second Comparison
        balcounts2 <- assay(ddsvbal2)
        balmeta2   <- colData(ddsvbal2)
        #Create Tables of Third Comparison
        balcounts3 <- assay(ddsvbal3)
        balmeta3   <- colData(ddsvbal3)
        #Create new DESEQ2 Objects
        ddsvbal1 <- DESeqDataSetFromMatrix(countData = balcounts1,
                                      colData = balmeta1,
                                      design= ~ three_groups)
        ddsvbal2 <- DESeqDataSetFromMatrix(countData = balcounts2,
                                      colData = balmeta2,
                                      design= ~ three_groups)
        ddsvbal3 <- DESeqDataSetFromMatrix(countData = balcounts3,
                                    colData = balmeta3,
                                    design= ~ three_groups)
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsvbal), 1, gm_mean)        
        #Estimate Factors of DESeq Object
        ddsvbal <- estimateSizeFactors(ddsvbal, geoMeans = geoMeans)
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsvbal1), 1, gm_mean)
        #Estimate Factors of DESeq Object
        ddsvbal1 <- estimateSizeFactors(ddsvbal1, geoMeans = geoMeans)
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsvbal2), 1, gm_mean)
        #Estimate Factors of DESeq Object
        ddsvbal2 <- estimateSizeFactors(ddsvbal2, geoMeans = geoMeans)
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsvbal3), 1, gm_mean)
        #Estimate Factors of DESeq Object
        ddsvbal3 <- estimateSizeFactors(ddsvbal3, geoMeans = geoMeans)
        #Variance Stablize the data
        vsdvbal <- varianceStabilizingTransformation(ddsvbal)
        #----------------------
        ##PCOA
        #----------------------
        #Create Distance Matrix
        vsdvbal0 <- ifelse(assay(vsdvbal)<0,0,assay(vsdvbal))
        vegdist   = vegdist(t(vsdvbal0), method="bray")
        #Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
        CmdScale <- cmdscale(vegdist, k =10)
        #calculated Sample variance for each PC
        vars <- apply(CmdScale, 2, var)
        #Create Variable with the Percent Variance
        percentVar <- round(100 * (vars/sum(vars)))
        #Merge PC Data with MetaData
        require(data.table)
        newResults <- merge(x = CmdScale, y = colData(vsdvbal), by = "row.names", all.x = TRUE)
        #Rename Variables for PC1 and PC2
        colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
        colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
        colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
        #Calculate the Centroid Value
        centroids <- aggregate(cbind(PC1,PC2)~ three_groups,data= newResults, mean)
        #Merge the Centroid Data into the PCOA Data
        newResults <- merge(newResults,centroids,by="three_groups",suffixes=c("",".centroid"))
        #Create Table for Statistics    
        data.adonis <- data.frame(colData(vsdvbal))
        #Run the Statistics
        samplepermanova <- adonis(vegdist ~ three_groups, data.adonis)
        samplepermanova <- as.data.frame(samplepermanova$aov.tab)
        samplepermanova <- samplepermanova$'Pr(>F)'[1]
        #PLOT IT
            ggsave(filename=paste0(counts,".BRAY_vsd_three_groups_PERMANOVA_",samplepermanova,".pdf"),
            ggplot(newResults, aes(PC1, PC2,color=three_groups)) + # Graph PC1 and PC2
            xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
            ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
            scale_color_manual(values=c("#D01C8B","#FFD479","#4DAC26" )) + 
            geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= three_groups))+ 
            geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Dead",">28 Days","<28 Days")), size=10) + #BKG DFW Lung.Tissue.In Lung.Tissue.UnIn Mock
            geom_point(data=newResults,aes(color=three_groups),size=5,alpha=0.5) + # Set the size of the points
            theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
            panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
            panel.grid.minor = element_blank(),strip.background=element_blank(),
            axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),
            axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),
            plot.margin=unit(c(1,1,1,1),"line"), legend.position="none"),
            height = 10, width = 10)
        #----------------------
        ##Alpha Diversity
        #----------------------
        #Calcultes Shannon Diversity
        ddsvbal$Shannon = diversity(vsdvbal0, index = "shannon", MARGIN = 2, base = exp(1))
        #Convert to data frame for ggplot
        shannon = as.data.frame(colData(ddsvbal))
        #Set Order Of Figure
        shannon$or <-ifelse(shannon$three_groups=="Less_Than_28_days_on_vent", 1,NA)
        shannon$or <-ifelse(shannon$three_groups=="Greater_Than_28_days_on_vent",2 ,shannon$or)
        shannon$or <-ifelse(shannon$three_groups=="Dead",3,shannon$or)
        #Make Sure Shannon is Numeric
        shannon$Shannon <- as.numeric(as.character(shannon$Shannon))
        #Check Statistics
        sampleshannon <- kruskal.test(Shannon ~ three_groups, data = shannon)
        sampleshannon <- sampleshannon$p.value
        #PLOT IT
            ggsave(filename=paste0(counts,".SHANNON_vsd_three_groups_Pvalue_",sampleshannon,".pdf"),
            ggplot(shannon, aes(x= reorder(three_groups, +or), y=Shannon, fill=three_groups)) + 
            stat_boxplot(geom ='errorbar', width=0.1)+
            geom_boxplot(outlier.shape = NA, width=0.5)+
            geom_jitter(shape=1, position=position_jitter(0.2))+
            scale_fill_manual(values=c("#D01C8B","#FFD479","#4DAC26")) + 
            scale_x_discrete(labels = c('<28 Days','>28 Days','Dead'))+ 
            ylab("Shannon Diversity") + 
            xlab("Three groups")+
            theme,
            height = 7, width = 5)
        #Check Statistics
        kruskal.test(Shannon ~ three_groups, data = shannon)
        #----------------------
        ##Differential Analysis
        #----------------------
        #DropLevels
        ddsvbal1$three_groups <- droplevels(ddsvbal1$three_groups)
        ddsvbal2$three_groups <- droplevels(ddsvbal2$three_groups)
        ddsvbal3$three_groups <- droplevels(ddsvbal3$three_groups)
        #Set Reference
        ddsvbal1$three_groups <- relevel(ddsvbal1$three_groups, ref ="Less_Than_28_days_on_vent")
        ddsvbal2$three_groups <- relevel(ddsvbal2$three_groups, ref ="Less_Than_28_days_on_vent")
        ddsvbal3$three_groups <- relevel(ddsvbal3$three_groups, ref ="Greater_Than_28_days_on_vent")
        #Run the differential Analysis: Lung Cancer Vs Wild Type --> positive is upregulated in Lung Cancer; Negative is down regulated
        ddsvbal1  <- DESeq(ddsvbal1)
        ddsvbal2  <- DESeq(ddsvbal2)
        ddsvbal3  <- DESeq(ddsvbal3)
        #Output Result Table
        res1     <- results(ddsvbal1, cooksCutoff=FALSE)
        res2     <- results(ddsvbal2, cooksCutoff=FALSE)
        res3     <- results(ddsvbal3, cooksCutoff=FALSE)
        #----------------------
        ##TABLES
        #----------------------
        #Get Assay Data For Compairson 1
        GenusData <-as.data.frame(assay(ddsvbal1)) #pruned to selected Genuses based on abundance
        #Create Relative Abundance Table
        df <-
            GenusData %>% 
                rownames_to_column('gs') %>%
                group_by(gs) %>% 
                summarise_all(funs(sum)) %>%
                mutate_if(is.numeric, funs(./sum(.))) %>%
                column_to_rownames('gs')
        #Get the ColData for Each Comparison
        coldata.1 <- coldata2[coldata2$three_groups=="Less_Than_28_days_on_vent",] %>%
                    select(Study_Linked_ID)
        coldata.2 <- coldata2[coldata2$three_groups=="Dead",] %>%
                    select(Study_Linked_ID)
        #keep Count data only for each comparison
        needed<-which(colnames(df) %in% rownames(coldata.1))    
        df.1 <- df[,needed]
        needed2<-which(colnames(df) %in% rownames(coldata.2))    
        df.2 <- df[,needed2]
        #Convert Resuts table into a data.frame
        res1 <- as.data.frame(res1)
        #decide what otu to save 
        otu.to.save <-as.character(rownames(res1))
        #from relative table we should get the mean across the row of the otu table
        df.1.meanRA <- rowMeans(df.1)
        df.2.meanRA <- rowMeans(df.2)
        #need to subset AND reorder just the otus that we have 
        df.1.meanRA.save <- df.1.meanRA[otu.to.save]
        df.2.meanRA.save <- df.1.meanRA[otu.to.save]
        #add the abundnace data for the res dataframe
        res1$abundance.1 <- df.1.meanRA.save
        res1$abundance.2 <- df.1.meanRA.save
        #Set Names of Results Table
        res1 <- setNames(cbind(rownames(res1), res1, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
        #Get Assay Data For Compairson 2
        GenusData <-as.data.frame(assay(ddsvbal2)) #pruned to selected Genuses based on abundance
        #Create Relative Abundance Table
        df <-
            GenusData %>% 
                rownames_to_column('gs') %>%
                group_by(gs) %>% 
                summarise_all(funs(sum)) %>%
                mutate_if(is.numeric, funs(./sum(.))) %>%
                column_to_rownames('gs')
        #Get the ColData for Each Comparison
        coldata.1 <- coldata2[coldata2$three_groups=="Less_Than_28_days_on_vent",] %>%
                    select(Study_Linked_ID)
        coldata.2 <- coldata2[coldata2$three_groups=="Greater_Than_28_days_on_vent",] %>%
                    select(Study_Linked_ID)
        #keep Count data only for each comparison
        needed<-which(colnames(df) %in% rownames(coldata.1))    
        df.1 <- df[,needed]
        needed2<-which(colnames(df) %in% rownames(coldata.2))    
        df.2 <- df[,needed2]
        #Convert Resuts table into a data.frame
        res2 <- as.data.frame(res2)
        #decide what otu to save 
        otu.to.save <-as.character(rownames(res2))
        #from relative table we should get the mean across the row of the otu table
        df.1.meanRA <- rowMeans(df.1)
        df.2.meanRA <- rowMeans(df.2)
        #need to subset AND reorder just the otus that we have 
        df.1.meanRA.save <- df.1.meanRA[otu.to.save]
        df.2.meanRA.save <- df.1.meanRA[otu.to.save]
        #add the abundnace data for the res dataframe
        res2$abundance.1 <- df.1.meanRA.save
        res2$abundance.2 <- df.1.meanRA.save
        #Set Names of Results Table
        res2 <- setNames(cbind(rownames(res2), res2, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
        #Get Assay Data For Compairson 3
        GenusData <-as.data.frame(assay(ddsvbal3)) #pruned to selected Genuses based on abundance
        #Create Relative Abundance Table
        df <-
            GenusData %>% 
                rownames_to_column('gs') %>%
                group_by(gs) %>% 
                summarise_all(funs(sum)) %>%
                mutate_if(is.numeric, funs(./sum(.))) %>%
                column_to_rownames('gs')
        #Get the ColData for Each Comparison
        coldata.1 <- coldata2[coldata2$three_groups=="Greater_Than_28_days_on_vent",] %>%
                    select(Study_Linked_ID)
        coldata.2 <- coldata2[coldata2$three_groups=="Dead",] %>%
                    select(Study_Linked_ID)
        #keep Count data only for each comparison
        needed<-which(colnames(df) %in% rownames(coldata.1))    
        df.1 <- df[,needed]
        needed2<-which(colnames(df) %in% rownames(coldata.2))    
        df.2 <- df[,needed2]
        #Convert Resuts table into a data.frame
        res3 <- as.data.frame(res3)
        #decide what otu to save 
        otu.to.save <-as.character(rownames(res3))
        #from relative table we should get the mean across the row of the otu table
        df.1.meanRA <- rowMeans(df.1)
        df.2.meanRA <- rowMeans(df.2)
        #need to subset AND reorder just the otus that we have 
        df.1.meanRA.save <- df.1.meanRA[otu.to.save]
        df.2.meanRA.save <- df.1.meanRA[otu.to.save]
        #add the abundnace data for the res dataframe
        res3$abundance.1 <- df.1.meanRA.save
        res3$abundance.2 <- df.1.meanRA.save
        #Set Names of Results Table
        res3 <- setNames(cbind(rownames(res3), res3, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
        #Write Tables of Differential Analysis
        write.table(res1,file=paste0(counts,".Less_Than_28D_vs_Dead_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        write.table(res2,file=paste0(counts,".Less_Than_28D_vs_Greater_Than_28D_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        write.table(res3,file=paste0(counts,".Greater_Than_28D_vs_Dead_all.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        #Remove Any Data without LOGFC data
        res1 <- res1[!is.na(res1$logFC),]
        res2 <- res2[!is.na(res2$logFC),]
        res3 <- res3[!is.na(res3$logFC),]
        #Match the three tables with the same Taxa
        needed<-which(res1$Gene.symbol %in% res3$Gene.symbol)    
        res1 <- res1[needed,]
        needed<-which(res2$Gene.symbol %in% res1$Gene.symbol)  
        res2 <- res2[needed,]
        needed<-which(res3$Gene.symbol %in% res1$Gene.symbol)  
        res3 <- res3[needed,]
        needed<-which(res1$Gene.symbol %in% res2$Gene.symbol)    
        res1 <- res1[needed,]
        needed<-which(res3$Gene.symbol %in% res2$Gene.symbol)    
        res3 <- res3[needed,]
        # Reorder Results based on FDR for comparison 1
        res1 = res1[order(res1$adj.P.Val, na.last = NA), ]
        # Keep only top 10
        res1order <- res1 %>% slice(1:10)
        # Reorder Results based on FDR for comparison 2
        res2 = res2[order(res2$adj.P.Val, na.last = NA), ]
        # Keep only top 10
        res2order <- res2 %>% slice(1:10)
        # Reorder Results based on FDR for comparison 3
        res3 = res3[order(res3$adj.P.Val, na.last = NA), ]
        # Keep only top 10
        res3order <- res3 %>% slice(1:10)
        #Combine order of three tables
        resorder <- rbind(res1order,res2order,res3order)
        #Remove any duplicates
        resorder <- resorder[!duplicated(resorder$Gene.symbol), ]
        #Keep only matching Taxa for all three tables
        res1 <- res1[res1$Gene.symbol %in% resorder$Gene.symbol,]
        #Set The order
        res1 <- res1[ order(match(res1$Gene.symbol, resorder$Gene.symbol)), ]
        #Create Variable for order
        res1 <- res1 %>% mutate(start = 1:n())
        #Keep only matching Taxa for all three tables
        res2 <- res2[res2$Gene.symbol %in% resorder$Gene.symbol,]
        #Set The order
        res2 <- res2[ order(match(res2$Gene.symbol, resorder$Gene.symbol)), ]
        #Keep only matching Taxa for all three tables
        res2 <- res2 %>% mutate(start = 1:n())
        #Keep only matching Taxa for all three tables
        res3 <- res3[res3$Gene.symbol %in% resorder$Gene.symbol,]
        #Set The order
        res3 <- res3[ order(match(res3$Gene.symbol, resorder$Gene.symbol)), ]
        #Keep only matching Taxa for all three tables
        res3 <- res3 %>% mutate(start = 1:n())
        #Set Variable for the three comparisson
        res1$group <- "A"
        res2$group <- "B"
        res3$group <- "C"
        #Convert Important columns to Numeric
        res1$adj.P.Val <-   as.numeric(as.character(res1$adj.P.Val))
        res1$logFC <-       as.numeric(as.character(res1$logFC))
        res1$abundance.1 <- as.numeric(as.character(res1$abundance.1))
        res1$abundance.2 <- as.numeric(as.character(res1$abundance.2))
        res2$adj.P.Val <-   as.numeric(as.character(res2$adj.P.Val))
        res2$logFC <-       as.numeric(as.character(res2$logFC))
        res2$abundance.1 <- as.numeric(as.character(res2$abundance.1))
        res2$abundance.2 <- as.numeric(as.character(res2$abundance.2))
        res3$adj.P.Val <-   as.numeric(as.character(res3$adj.P.Val))
        res3$logFC <-       as.numeric(as.character(res3$logFC))
        res3$abundance.1 <- as.numeric(as.character(res3$abundance.1))
        res3$abundance.2 <- as.numeric(as.character(res3$abundance.2))
        #Bind the three Comparison TAbles
        resy <- rbind(res1,res2,res3)
        #Create Variable for Color based on Comparison, FDR and LOGFC
        resy$col <- ifelse(resy$group=="A" & resy$adj.P.Val<0.2 & resy$logFC>0, "A",
                    ifelse(resy$group=="A" & resy$adj.P.Val<0.2 & resy$logFC<0, "B",
                    ifelse(resy$group=="B" & resy$adj.P.Val<0.2 & resy$logFC>0, "C",
                    ifelse(resy$group=="B" & resy$adj.P.Val<0.2 & resy$logFC<0, "B",
                    ifelse(resy$group=="C" & resy$adj.P.Val<0.2 & resy$logFC>0, "A",
                    ifelse(resy$group=="C" & resy$adj.P.Val<0.2 & resy$logFC<0, "C","D"))))))
        #PLOT IT
            ggsave(filename=paste0(counts,".BAL_three_groups_DESEQ2.pdf"),
            ggplot(resy, aes(y=reorder(Gene.symbol,-start), x=logFC,fill=col)) +
            facet_grid(~ group, scales = "free_y")+
            geom_point(size = ifelse(resy$adj.P.Val<0.2 & resy$logFC>0, 2000 * resy$abundance.2, ifelse(resy$adj.P.Val<0.2 & resy$logFC<0, 2000 * resy$abundance.1,5)),color="black",alpha=0.8,shape=21)+
            scale_fill_manual(values=c("#D01C8B","#4DAC26","#FFD479","white"))+
            scale_size_continuous(range=c(1, 27),guide=FALSE)+
            theme(panel.background = element_blank(),
                panel.border=element_rect(fill=NA),
                panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
                panel.grid.minor = element_blank(),
                strip.background=element_blank(),
                axis.title=element_text(size=20,face="bold"),
                axis.text.x=element_text(colour="black", size=18, face="bold"),
                axis.text.y=element_text(colour="black",face="bold",size=10),
                axis.ticks=element_line(colour="black"),
        		legend.background = element_rect(color=NA))+
            xlab("") +
            ylab("")+
            xlim(-2,7)+
        	geom_vline(xintercept=0, color="red",linetype="dashed")+
        	guides(fill=FALSE),
            width=20, height=5)

}

