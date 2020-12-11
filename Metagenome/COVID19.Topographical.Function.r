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
deseq <- function(metadata,counts,oldcounts,sampletype,comparison) {
        #Load MetaData
        coldata <- read.delim2(metadata, sep="\t")
        rownames(coldata) <- coldata$GTC_Sequence_ID
        #Fix Three Groups
        coldata$three_groups <- ifelse(coldata$X.3_groups=="\xb2.28.day.on.vent","Less_Than_28_days_on_vent",
                                ifelse(coldata$X.3_groups==">28.day.on.vent","Greater_Than_28_days_on_vent",
                                as.character(coldata$X.3_groups)))
        #load Count Data
        mycounts <-read.delim2(paste0(counts,".txt"), sep="\t", row.names=1)
        mycountsold <-read.delim2(paste0(oldcounts,".txt"), sep="\t", row.names=1)
        #Find Unique Label for SampleID (multiple LANEs)
        prefixes = unique(unlist(lapply(colnames(mycounts), function(x) { strsplit(x, '_')[[1]][[1]]})))
        #Creat Function to do sums of same IDs
        f = function(prefix) {
          return (rowSums(subset(mycounts, select=grep(prefix, colnames(mycounts)))))
        }
        #Run Function
        m = matrix(unlist(lapply(prefixes, f)), nrow=nrow(mycounts))
        #Fix the colnames and Rownames
        colnames(m) = prefixes
        rownames(m) = rownames(mycounts)
        #Convert BAck to DataFrame
        mycounts <- as.data.frame(m)
        #Fix the Names to match the Coldata
        names(mycounts) <- gsub(x = names(mycounts), pattern = "\\.", replacement = "_")
        rownames(mycounts) <- gsub(x = rownames(mycounts), pattern = " ", replacement = "_")  
        rownames(mycountsold) <- gsub(x = rownames(mycountsold), pattern = " ", replacement = "_")  
        #Find matching sample ID for both datasets
        needed<-which(rownames(coldata) %in% colnames(mycounts))    
        neededold<-which(rownames(coldata) %in% colnames(mycountsold))    
        #keep only matching IDs from count data
        coldata2<-coldata[needed,]
        coldataold<-coldata[neededold,]
        #Order Meta Data by SampleId
        coldata2 <- coldata2[order(coldata2$GTC_Sequence_ID),]
        coldataold <- coldataold[order(coldataold$GTC_Sequence_ID),]
        #keep only matching IDs from count data
        wanted<-which(colnames(mycounts) %in% rownames(coldata))    
        wantedold<-which(colnames(mycountsold) %in% rownames(coldataold))    
        #keep only matching IDs from count data
        mycounts2<-mycounts[,wanted]
        mycounts2old<-mycountsold[,wantedold]
        #Order Count Data by SampleID
        mycounts2 <-mycounts2[, order(colnames(mycounts2))]
        mycounts2old <-mycounts2old[, order(colnames(mycounts2old))]
        #Convert any NAs to 0
        mycounts2[is.na(mycounts2)] <- 0
        mycounts2old[is.na(mycounts2old)] <- 0
        #Copy of Count Table
        mycounts3 <- mycounts2
        mycounts3old <- mycounts2old
        #Convert Count Table into a Numeic Data Frame
        d1 = data.frame(lapply(mycounts3, function(x) as.numeric(as.character(x))),
                           check.names=F, row.names = rownames(mycounts3))
        d2 = data.frame(lapply(mycounts3old, function(x) as.numeric(as.character(x))),
                           check.names=F, row.names = rownames(mycounts3old))
        #Convert Data to Integers to Run DESEq
        d1[] <- lapply(d1, as.integer)
        d2[] <- lapply(d2, as.integer)
        #Combine Meta and Count Data
        ddsvbkg <- DESeqDataSetFromMatrix(countData = d1,
                              colData = coldata2,
                              design= ~ Sample.Type)
        #keep on BKG
        ddsvbkg <- ddsvbkg[, ddsvbkg$Sample.Type %in% c("BKG")]
        #Get seperate Tables
        bkgcountsnew <- assay(ddsvbkg)
        bkgmetanew   <- colData(ddsvbkg)
        #Fix Rownames
        rownames(bkgcountsnew) <- gsub(x = rownames(bkgcountsnew), pattern = " ", replacement = "_")  
        #Combine Old Meta and Count Data
        ddsvbalua <- DESeqDataSetFromMatrix(countData = d2,
                              colData = coldataold,
                              design= ~ Sample.Type)   
        #keep on BKG
        ddsvbalua <- ddsvbalua[, ddsvbalua$Sample.Type %in% c("BAL","UA")]
        #Get seperate Tables     
        baluacounts <- assay(ddsvbalua)
        baluameta   <- colData(ddsvbalua)
        #balmeta$batch <- "Old"
        needed<-which(rownames(baluacounts) %in% rownames(bkgcountsnew))    
        baluacounts<-baluacounts[needed,]
        needed<-which(rownames(bkgcountsnew) %in% rownames(baluacounts))    
        bkgcountsnew<-bkgcountsnew[needed,]
        #Merge DataSets
        combinedcount <- cbind(baluacounts,bkgcountsnew)
        combinedmeta  <- rbind(baluameta,bkgmetanew)
        #Fix one variable without Group
        combinedmeta$three_groups <-ifelse(is.na(combinedmeta$three_groups),"Greater_Than_28_days_on_vent",as.character(combinedmeta$three_groups))
        #Create DESEQ2 Pbkect
        ddsv <- DESeqDataSetFromMatrix(countData = combinedcount,
                            colData = combinedmeta,
                            design= ~ three_groups)
        write.table(assay(ddsv),file=paste0(counts,".Merged_Count_Table.txt"), sep="\t", col.names = NA, row.names = TRUE)
        write.table(colData(ddsv),file=paste0(counts,".Merged_MetaData.txt"), sep="\t", col.names = NA, row.names = TRUE)
        #remove 0 coutns  
        idx <- which.max(colSums(counts(ddsv) == 0))
        ddsv <- ddsv[ , -idx]
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
        ##Differential Analysis
        #----------------------
        #Create Deseq object for BALUA analysis
        ddsvbal <- DESeqDataSetFromMatrix(countData = combinedcount,
                              colData = combinedmeta,
                              design= ~ Sample.Type)
        #Subset just the BAL
        ddsvbal <- ddsvbal[, ddsvbal$Sample.Type %in% c("BAL","UA")]
        #Covert Variable to Factor
        ddsvbal$Sample.Type <- as.factor(ddsvbal$Sample.Type)
        #Calculate geometric means prior to estimate size factor
        gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
        geoMeans = apply(counts(ddsvbal), 1, gm_mean)        
        #Estimate Factors of DESeq Object
        ddsvbal <- estimateSizeFactors(ddsvbal, geoMeans = geoMeans)
        #Variance Stablize the data
        vsdvbal <- varianceStabilizingTransformation(ddsvbal)
        #DropLevels
        ddsvbal$Sample.Type <- droplevels(ddsvbal$Sample.Type)
        #Set Reference
        ddsvbal$Sample.Type <- relevel(ddsvbal$Sample.Type, ref ="UA")
        #Run the differential Analysis: Lung Cancer Vs Wild Type --> positive is upregulated in Lung Cancer; Negative is down regulated
        ddsvbal  <- DESeq(ddsvbal)
        #Output Result Table
        res1     <- results(ddsvbal, cooksCutoff=FALSE)
        #----------------------
        ##TABLES
        #----------------------
        #Get Assay Data For Compairson 1
        GenusData <-as.data.frame(assay(ddsvbal)) #pruned to selected Genuses based on abundance
        #Create Relative Abundance Table
        df <-
            GenusData %>% 
                rownames_to_column('gs') %>%
                group_by(gs) %>% 
                summarise_all(funs(sum)) %>%
                mutate_if(is.numeric, funs(./sum(.))) %>%
                column_to_rownames('gs')
        #Get the ColData for Each Comparison
        coldata.1 <- coldata2[coldata2$Sample.Type=="UA",] %>%
                    select(GTC_Sequence_ID)
        coldata.2 <- coldata2[coldata2$Sample.Type=="BAL",] %>%
                    select(GTC_Sequence_ID)
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
        df.2.meanRA.save <- df.2.meanRA[otu.to.save]
        #add the abundnace data for the res dataframe
        res1$abundance.1 <- df.1.meanRA.save
        res1$abundance.2 <- df.2.meanRA.save
        #Set Names of Results Table
        res1 <- setNames(cbind(rownames(res1), res1, row.names = NULL), c("Gene.symbol","baseMean", "logFC", "lfcSE", "stat", "pvalue", "adj.P.Val","abundance.1","abundance.2")) 
        #Get Assay Data For Compairson 2
        #Write Tables of Differential Analysis
        write.table(res1,file=paste0(counts,".UA_vs_BAL.txt"), sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
        #Remove Any Data without LOGFC data
        #=========================================================
        #------------VOLCANO PLOT
        #=========================================================
        # Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
        res1$sig <- -log10(res1$adj.P.Val)
        sum(is.infinite(res1$sig))
        res1[is.infinite(res1$sig),"sig"] <- 350
        ## Set Color Gradient
        cols <- densCols(res1$logFC, res1$sig)
        cols[res1$pvalue ==0] <- "purple"
        cols[res1$logFC > 0 & res1$adj.P.Val < alpha ] <- "red"
        cols[res1$logFC < 0 & res1$adj.P.Val < alpha ] <- "green"
        res1$pch <- 19
        res1$pch[res1$pvalue ==0] <- 6
            ggsave(filename=paste0(counts,".DESEQ2_UA_VS_BAL_FDR0.05.pdf"),
            ggplot(res1, aes(x = logFC, y = sig,label=Gene.symbol)) +
            geom_point(color=cols, size = ifelse(res1$logFC>=1 & res1$adj.P.Val < alpha, 200 * res1$abundance.2, ifelse(res1$logFC<=-1 & res1$adj.P.Val < alpha, 200 * res1$abundance.1,1)),alpha=0.6) + #Chose Colors for dots
            geom_text_repel(aes(label=ifelse(res1$logFC<(-3) & res1$adj.P.Val < alpha , as.character(res1$Gene.symbol),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
            geom_text_repel(aes(label=ifelse(res1$logFC>3 & res1$adj.P.Val < alpha , as.character(res1$Gene.symbol),'')),size=3,force=25,segment.colour="grey",segment.alpha=0.5) +
            theme(legend.position = "none") +
            geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
            xlab("Effect size: log2(fold-change)") +
            ylab("-log10(adjusted p-value)") + 
            #ylim(0,20)+
            theme, 
            width=5, height=5)
}

