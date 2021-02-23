#Figures from the Paper
#Figure 3B
ssh -Y sulaii01@bigpurple.nyumc.org

#Connect to the Node
srun -p cpu_medium --nodes=1 --tasks-per-node=2 --cpus-per-task=1 --mem-per-cpu=40G -t 04:00:00  --pty bash

#Load Modules
module load r/3.5.1
module load gsl/1.15

#Open R
R

#Load Libraries
library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(eulerr)
library(tibble)

#Set Theme for Figures
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")

##################################################3
#SPECIES
##################################################

#Less than 28 Days versus Dead
#Load Metatranscriptome Data
setwd("/gpfs/data/segallab/COVID19/MultiOmics/Bacterial/")
x <- read.table("metatrans_bacteria.Less_Than_28D_vs_Dead_all.txt", header=T, sep="\t", row.names=1)
y <- read.table("metatrans_bacteria.Less_Than_28D_vs_Greater_Than_28D_all.txt", header=T, sep="\t", row.names=1)
z <- read.table("metatrans_bacteria.Greater_Than_28D_vs_Dead_all.txt", header=T, sep="\t", row.names=1)

#Load Metagenome Data
#setwd("/gpfs/data/segallab/COVID19/Metagenome/Bacterial/Species/")
a <- read.table("metagenome_bacteria.Less_Than_28D_vs_Dead_all.txt", header=T, sep="\t", row.names=1)
b <- read.table("metagenome_bacteria.Less_Than_28D_vs_Greater_Than_28D_all.txt", header=T, sep="\t", row.names=1)
c <- read.table("metagenome_bacteria.Greater_Than_28D_vs_Dead_all.txt", header=T, sep="\t", row.names=1)

#//Metagenome//
#Remove NAs
a$Genus <- a$Gene.symbol
a$Genus <- gsub("_", " ",a$Genus)
a$Genus <- gsub("_", " ",a$Genus)

#a <- a %>% 
#  rowwise() %>% 
#  mutate(abundance = median(c(abundance.1, abundance.2), na.rm = TRUE))
data.a <- a[!is.na(a$Genus),]
data.a <- data.a[!is.na(data.a$logFC),]
data.a <- data.a[!is.na(data.a$adj.P.Val),]
data.a$seq <- "Metagenome"
data.a$compare <- "A"
data.a$sig <- -log10(data.a$adj.P.Val)

b$Genus <- b$Gene.symbol
b$Genus <- gsub("_", " ",b$Genus)
b$Genus <- gsub("_", " ",b$Genus)

data.b <- b[!is.na(b$Genus),]
data.b <- data.b[!is.na(data.b$logFC),]
data.b <- data.b[!is.na(data.b$adj.P.Val),]
data.b$seq <- "Metagenome"
data.b$compare <- "B"
data.b$sig <- -log10(data.b$adj.P.Val)

c$Genus <- c$Gene.symbol
c$Genus <- gsub("_", " ",c$Genus)
c$Genus <- gsub("_", " ",c$Genus)

data.c <- c[!is.na(c$Genus),]
data.c <- data.c[!is.na(data.c$logFC),]
data.c <- data.c[!is.na(data.c$adj.P.Val),]
data.c$seq <- "Metagenome"
data.c$compare <- "C"
data.c$sig <- -log10(data.c$adj.P.Val)

#Label Sequencing Type
#data.a$seq <- "Metagenome"
##convert adjusted pvalue to log10
#data.a$sig <- -log10(data.a$adj.P.Val)
##Calculate Median IQR and N of Padg
#data.a <- setDT(data.a)[,list(Number=median(abundance),medianfc=median(logFC), median=as.numeric(median(adj.P.Val)), iqr=as.numeric(quantile(logFC, probs=.75)),iqr1=as.numeric(quantile(logFC, probs=.25))), by=c("Genus")]
##Label Sequencing Type
#data.a$seq <- "Metagenome"
##convert adjusted pvalue to log10
#data.a$sig <- -log10(data.a$median)
##Make sure no values become infinity
#sum(is.infinite(data.a$sig))

#//Metatranscriptome//
#Remove NAs
x$Genus <- x$Gene.symbol
#x <- x %>% 
#  rowwise() %>% 
#  mutate(abundance = median(c(abundance.1, abundance.2), na.rm = TRUE))
data.x <- x[!is.na(x$Genus),]
data.x <- data.x[!is.na(data.x$logFC),]
data.x <- data.x[!is.na(data.x$adj.P.Val),]
data.x$seq <- "Metatranscriptome"
data.x$compare <- "A"
data.x$sig <- -log10(data.x$adj.P.Val)

y$Genus <- y$Gene.symbol
#x <- x %>% 
#  rowwise() %>% 
#  mutate(abundance = median(c(abundance.1, abundance.2), na.rm = TRUE))
data.y <- y[!is.na(y$Genus),]
data.y <- data.y[!is.na(data.y$logFC),]
data.y <- data.y[!is.na(data.y$adj.P.Val),]
data.y$seq <- "Metatranscriptome"
data.y$compare <- "B"
data.y$sig <- -log10(data.y$adj.P.Val)

z$Genus <- z$Gene.symbol
#x <- x %>% 
#  rowwise() %>% 
#  mutate(abundance = median(c(abundance.1, abundance.2), na.rm = TRUE))
data.z <- z[!is.na(z$Genus),]
data.z <- data.z[!is.na(data.z$logFC),]
data.z <- data.z[!is.na(data.z$adj.P.Val),]
data.z$seq <- "Metatranscriptome"
data.z$compare <- "C"
#convert adjusted pvalue to log10
data.z$sig <- -log10(data.z$adj.P.Val)
#
##Calculate Median IQR and N of Padg
#data.x <- setDT(data.x)[,list(Number=median(abundance),medianfc=median(logFC), median=as.numeric(median(adj.P.Val)), iqr=as.numeric(quantile(logFC, probs=.75)),iqr1=as.numeric(quantile(adj.P.Val, probs=.25))), by=c("Genus")]
##Label Sequencing Type
#data.x$seq <- "Metatranscriptome"
##convert adjusted pvalue to log10
#data.x$sig <- -log10(data.x$median)
##Make sure no values become infinity
#sum(is.infinite(data.x$sig))

#Create Ranks
head(x)
res1 <- x %>% 
   dplyr::select(Genus, stat) %>% 
   na.omit() %>% 
   distinct() %>% 
   group_by(Genus) %>% 
   summarize(stat=mean(as.numeric(stat)))
res1

#Deframe Ranks
ranks <- deframe(res1)
head(ranks, 20)

#make sure numeric variables are numeric
a$logFC <- as.numeric(a$logFC)
a$adj.P.Val <- as.numeric(a$adj.P.Val)
a <- a[!is.na(a$adj.P.Val),]
#a <- a[a$padj <= 0.1,] # set pdj threshold
a$fcSign <- sign(a$logFC)
gmt.file<-c()
gmt.file$GenesetUp<-a[a$fcSign==1,]$Genus
gmt.file$GenesetDown<-a[a$fcSign==-1,]$Genus

#Run GSEA Analysis comparing Metatranscriptome to Metagenome
fgseaRes <- fgsea(pathways=gmt.file, stats=ranks, nperm=1000)

#Covert Output to a Data.Frame
gsea <- as.data.frame(fgseaRes)

#Spit the KOs into their own row
library(splitstackshape)
gsea <- cSplit(as.data.table(gsea)[, leadingEdge := gsub("[][\"]", "", leadingEdge)], 
       "leadingEdge", ",", "long")
#Remove the c(
gsea[] <- lapply(gsea, gsub, pattern='c\\(', replacement='')
#Remove the )
gsea[] <- lapply(gsea, gsub, pattern=')', replacement='')

#Count Number of Overlapping KOs
length(unique(gsea[["leadingEdge"]])) #135
#Count Number of KOs in Metatranscriptome Data
length(unique(x[["Genus"]])) #2572
#Count Number of KOs in Metagenome Data
length(unique(a[["Genus"]])) #904

#Metatranscriptome (58) is Number minus the overalp and Metagenome (16) is Number minus the overlap
fit <- euler(c(Metatranscriptome = 2437, Metagenome = 769, "Metatranscriptome&Metagenome" = 135)) 

#---------------
#----Figure 3A
#---------------
pdf("Metatranscriptome_vs_Metagenome_Less_than_28_Days_vs_Dead_SPECIES_VENN.pdf", height = 20, width = 25)
plot(fit, quantities=list(cex = 3.5),fills = list(fill = c("#7A81FF","#FFD479"), alpha = 0.5),
     labels = list(col = "black", cex=3.5))
dev.off()


#B
head(y)
res2 <- y %>% 
   dplyr::select(Genus, stat) %>% 
   na.omit() %>% 
   distinct() %>% 
   group_by(Genus) %>% 
   summarize(stat=mean(as.numeric(stat)))
res2

#Deframe Ranks
ranks2 <- deframe(res2)
head(ranks2, 20)

#make sure numeric variables are numeric
b$logFC <- as.numeric(b$logFC)
b$adj.P.Val <- as.numeric(b$adj.P.Val)
b <- b[!is.na(b$adj.P.Val),]
#a <- a[a$padj <= 0.1,] # set pdj threshold
b$fcSign <- sign(b$logFC)
gmt.file2<-c()
gmt.file2$GenesetUp<-  b[b$fcSign==1,]$Genus
gmt.file2$GenesetDown<-b[b$fcSign==-1,]$Genus

#Run GSEA Analysis comparing Metatranscriptome to Metagenome
fgseaRes2 <- fgsea(pathways=gmt.file2, stats=ranks2, nperm=1000)

#Covert Output to a Data.Frame
gsea2 <- as.data.frame(fgseaRes2)

#Spit the KOs into their own row
library(splitstackshape)
gsea2 <- cSplit(as.data.table(gsea2)[, leadingEdge := gsub("[][\"]", "", leadingEdge)], 
       "leadingEdge", ",", "long")
#Remove the c(
gsea2[] <- lapply(gsea2, gsub, pattern='c\\(', replacement='')
#Remove the )
gsea2[] <- lapply(gsea2, gsub, pattern=')', replacement='')


#Count Number of Overlapping KOs
length(unique(gsea2[["leadingEdge"]])) #533
#Count Number of KOs in Metatranscriptome Data
length(unique(y[["Genus"]])) #2572
#Count Number of KOs in Metagenome Data
length(unique(b[["Genus"]])) #3639

#Metatranscriptome (58) is Number minus the overalp and Metagenome (16) is Number minus the overlap
fit <- euler(c(Metatranscriptome = 2039, Metagenome = 3106, "Metatranscriptome&Metagenome" = 533)) 

#---------------
#----Figure 3A
#---------------
pdf("Metatranscriptome_vs_Metagenome_Greater_than_28_Days_vs_Less_than_28_Days_SPECIES_VENN.pdf", height = 20, width = 25)
plot(fit, quantities=list(cex = 3.5),fills = list(fill = c("#7A81FF","#FFD479"), alpha = 0.5),
     labels = list(col = "black", cex=3.5))
dev.off()


#C
head(z)
res3 <- z %>% 
   dplyr::select(Genus, stat) %>% 
   na.omit() %>% 
   distinct() %>% 
   group_by(Genus) %>% 
   summarize(stat=mean(as.numeric(stat)))
res3

#Deframe Ranks
ranks3 <- deframe(res3)
head(ranks3, 20)

#make sure numeric variables are numeric
c$logFC <- as.numeric(c$logFC)
c$adj.P.Val <- as.numeric(c$adj.P.Val)
c <- c[!is.na(c$adj.P.Val),]
c$fcSign <- sign(c$logFC)
gmt.file3<-c()
gmt.file3$GenesetUp<-  c[c$fcSign==1,]$Genus
gmt.file3$GenesetDown<-c[c$fcSign==-1,]$Genus

#Run GSEA Analysis comparing Metatranscriptome to Metagenome
fgseaRes3 <- fgsea(pathways=gmt.file3, stats=ranks3, nperm=1000)

#Covert Output to a Data.Frame
gsea3 <- as.data.frame(fgseaRes3)

#Spit the KOs into their own row
library(splitstackshape)
gsea3 <- cSplit(as.data.table(gsea3)[, leadingEdge := gsub("[][\"]", "", leadingEdge)], 
       "leadingEdge", ",", "long")
#Remove the c(
gsea3[] <- lapply(gsea3, gsub, pattern='c\\(', replacement='')
#Remove the )
gsea3[] <- lapply(gsea3, gsub, pattern=')', replacement='')



#Count Number of Overlapping KOs
length(unique(gsea3[["leadingEdge"]])) #371
#Count Number of KOs in Metatranscriptome Data
length(unique(z[["Genus"]])) #2572
#Count Number of KOs in Metagenome Data
length(unique(c[["Genus"]])) #3207

#Metatranscriptome (58) is Number minus the overalp and Metagenome (16) is Number minus the overlap
fit <- euler(c(Metatranscriptome = 2201, Metagenome = 2836, "Metatranscriptome&Metagenome" = 371)) 

#---------------
#----Figure 3A
#---------------
pdf("Metatranscriptome_vs_Metagenome_Greater_than_28_Days_vs_Dead_SPECIES_VENN.pdf", height = 20, width = 25)
plot(fit, quantities=list(cex = 3.5),fills = list(fill = c("#7A81FF","#FFD479"), alpha = 0.5),
     labels = list(col = "black", cex=3.5))
dev.off()

#keep only shared KO_Subclass
keepKOs <- gsea$leadingEdge
data.a <- data.a[data.a$Genus %in% keepKOs,]
data.x <- data.x[data.x$Genus %in% keepKOs,]
data.a <- data.a[data.a$Genus %in% data.x$Genus,]

keepKOs2 <- gsea2$leadingEdge
data.b <- data.b[data.b$Genus %in% keepKOs2,]
data.y <- data.y[data.y$Genus %in% keepKOs2,]
data.b <- data.b[data.b$Genus %in% data.y$Genus,]
data.y <- data.y[data.y$Genus %in% data.b$Genus,]

keepKOs3 <- gsea2$leadingEdge
data.c <- data.c[data.c$Genus %in% keepKOs3,]
data.z <- data.z[data.z$Genus %in% keepKOs3,]
data.z <- data.z[data.z$Genus %in% data.c$Genus,]
data.c <- data.c[data.c$Genus %in% data.z$Genus,]


save.image(file="Combined.OMICS2.RData")


needed<-which(data.a$Genus %in% data.b$Genus)    
data.a <- data.a[needed,]
needed<-which(data.b$Genus %in% data.a$Genus)  
data.b <- data.b[needed,]
needed<-which(data.c$Genus %in% data.b$Genus)  
data.c <- data.c[needed,]
needed<-which(data.a$Genus %in% data.c$Genus)  
data.a <- data.a[needed,]
needed<-which(data.b$Genus %in% data.c$Genus)  
data.b <- data.b[needed,]

data.x <- data.x[data.x$Genus %in% data.a$Genus,]
data.y <- data.y[data.y$Genus %in% data.b$Genus,]
data.z <- data.z[data.z$Genus %in% data.c$Genus,]

#data.x <- data.x %>% select(-abundance)
#data.b <- data.b %>% select(-fcSign)
#data.a <- data.a %>% select(-abundance)
#
#data.b <- data.b %>% select(-sig)
#data.c <- data.c %>% select(-sig)
#data.y <- data.y %>% select(-sig)
#data.z <- data.z %>% select(-sig)


datas <- rbind(data.a,data.b,data.c,data.x,data.y,data.z)
sigs <-  datas %>% filter(adj.P.Val<0.2) %>% distinct(Genus)
sigs <-  datas %>% filter(adj.P.Val<0.05) %>% distinct(Genus)

datass <- datas[datas$Genus %in% sigs$Genus,]
abundance <-  datass %>% filter(abundance.1>0.0005|abundance.2>0.0005) %>% distinct(Genus)

df <- datas[datas$Genus %in% abundance$Genus,]
df <- datas[datas$Genus %in% sigs$Genus,]

#merge 3 datasets
#data <- rbind(data.a,data.x)
df$medianfc <- df$logFC
df$median <- df$adj.P.Val

#change Order of dataframe to be order of 16S data by Significance
df <- df[with(df, order(seq, +sig, Genus)),]
#change Order of dataframe to be order of 16S data by logFC
df <- df[with(df, order(seq, +medianfc, Genus)),]

#Keep the Uniqe Pathways
df$Genus <- factor(df$Genus, levels = unique(df$Genus))

#Set the shape of the Figures
df$shape<- as.character(as.numeric(ifelse(df$seq=="Metatranscriptome",21,23)))

#Covert LOGFC to Log10
df$log10 <- df$medianfc*log10(2)
#Fix the value that is less than -4 and greater than 5
df$log10 <- ifelse(df$log10< -6, -2, df$log10)
df$log10 <- ifelse(df$log10>5, df$log10-1, df$log10)


#set variables for Colors
#df$col <- ifelse(df$seq=="Metagenome" & df$median<0.05, "A",
#            ifelse(df$seq=="Metatranscriptome" & df$median<0.05, "B", "C"))

#Subset Only Significant Genus Names
sigs <-  df %>% filter(median<0.05) %>% distinct(Genus)
df2 <- df[df$Genus %in% sigs$Genus,]

df$col <-   ifelse(df$compare=="A" & df$adj.P.Val<0.2 & df$logFC>0, "A",
            ifelse(df$compare=="A" & df$adj.P.Val<0.2 & df$logFC<0, "B",
            ifelse(df$compare=="B" & df$adj.P.Val<0.2 & df$logFC>0, "C",
            ifelse(df$compare=="B" & df$adj.P.Val<0.2 & df$logFC<0, "B",
            ifelse(df$compare=="C" & df$adj.P.Val<0.2 & df$logFC>0, "A",

#Make a Copy of dataframe   
df3 <- df

#Plot Figure
pdf("COVID19_BAL_Metatranscriptome_vs_Metagenome_SPECIES_Bubble_Chart_sig_all_V2.pdf", height = 15, width = 20)
    ggplot(df, aes(x=logFC, y=reorder(Genus,+logFC), fill=col,shape=seq)) +
    facet_grid(~ compare, scales = "free_y")+
	scale_shape_manual(values=c(23,21), guide=FALSE) +
    geom_point(size = ifelse(df$adj.P.Val<0.2 & df$logFC>0, 2000 * df$abundance.2, ifelse(df$adj.P.Val<0.2 & df$logFC<0, 2000 * df$abundance.1,5)),color="black",alpha=0.8)+
    #geom_point(aes(fill=col),color="black",alpha=0.8)+
	#geom_point(color="black")+
    scale_size_continuous(range=c(1, 27),guide=FALSE)+
    scale_fill_manual(values=c("#D01C8B","#4DAC26","#FFD479","white"))+
    #scale_x_reverse()+ #Flip x axis so it goes from least significant to most
    #scale_colour_gradient(low="blue", high="black")+ #Set Color for gradient
    theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_text(colour="black",face="bold",size=10),
        axis.ticks=element_line(colour="black"),
        #plot.margin=unit(c(1,1,1,1),"line"),
        #legend.position = c(.95, .05),
        #legend.justification = c("right", "bottom"),
		#legend.key = element_rect(colour = NA),
		legend.background = element_rect(color=NA))+
    #facet_grid(bb ~ .)+
    #scale_colour_continuous(guide = FALSE)
	#scale_fill_manual(values=c("#FFD479","#7A81FF","white"), guide=FALSE) +
	#scale_color_manual(values=c( "#FF2F92","#FFD479","#7A81FF","white")) +
    #scale_x_continuous(limits=c(0,NA),expand=c(0,0),breaks=c(0,45:60), labels=c(0,45:60))+
    #annotate(geom = "segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)+
    xlab("") +
    ylab("")+
    #labs( size = "Number of KOs", color="Median Log P Value" ) +
	geom_vline(xintercept=0, color="red",linetype="dashed")+
	guides(fill=FALSE)
	#theme_bw()
dev.off()

