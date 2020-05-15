#' ---
#' title: "Specialty Figures"
#' output: word_document
#' ---
#' 
## ----setup, cache=FALSE, include=FALSE, warning = FALSE------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
opts_chunk$set(fig.width=4, fig.height=4, dpi=300, warning=FALSE, message=FALSE)

## ----include = FALSE-----------------------------------------------------
oldPar = par()

#' 
#' ## I. Menu
#' 
#' 
#' 
#' ##### --------------------------------------------------------------------------------------------------------
#' 
#' ### 1. Set Environment  
#' 
## ----warning = FALSE-----------------------------------------------------
#setwd("/Users/bswhemployee/Uthra_work/Rcourse_NIH/plots/")

# set working directory 

#? rm
rm(list = ls())


load("helperFunctions.RData")
load('NeveG.RData')
Data <- NeveG_data[, -(1:4)]
Info <- NeveG_data[,1:4]

str(NeveG_data) #50*13

#' 
#' ### 2. Word cloud
#' 
#' A word cloud is a text mining method that highlights the most frequently used words in a text. 
#' 
## ------------------------------------------------------------------------

library(wordcloud); 
library(tm); 
library("SnowballC"); 
library(RColorBrewer)

#SnowballC: implements Porter's word stemming algorithm 
# for collapsing words to a common root to aid comparison

# get a text file
fileInput <- "wordcloud_input.txt"
input <- readLines(fileInput)

# format as a corpus: 
#collections of documents containing (natural language) text

docs <- Corpus(VectorSource(input))
#inspect(docs)

# text transformations
toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))

#use 'tospace' functon above to replace special characters (/, @, |) into blank space 
# there are some warining messages bu that is OK. nothing is dropped or changed 

docs <- tm_map(docs, toSpace, "/")
docs <- tm_map(docs, toSpace, "@")
docs <- tm_map(docs, toSpace, "\\|")

# more transformation with  tm_map() function 
# Convert the text to lower case, so "data"= "Data"
docs <- tm_map(docs, content_transformer(tolower))

# Remove numbers
docs <- tm_map(docs, removeNumbers)

# Remove english common stopwords # such i, me, as, be # find a list online
docs <- tm_map(docs, removeWords, stopwords("english"))

# Remove your own stop word
# specify your own stopwords as a character vector
docs <- tm_map(docs, removeWords, c("ttt", "aaaa"))   


# Remove punctuations
docs <- tm_map(docs, removePunctuation)

# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)

# Text stemming, stem words from the common root. 
# Stem words in a text document using Porter's stemming algorithm
# This means that all the words are converted to their stem (Ex: learning -> learn, walked -> walk, etc.). 
#
#docs <- tm_map(docs, stemDocument)

# Create a term-document matrix
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)

v <- sort(rowSums(m), decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)

# get the words and frequency
# head(d, 10)

set.seed(123)
wordcloud(words = d$word, freq = d$freq, min.freq = 2,
          max.words=200, random.order=FALSE, rot.per=0.2, 
          colors=brewer.pal(8, "Dark2")) #Accent

# try again for different random layout

## Check the run a ML wiki text to make another World Cloud ##

input <- readLines(file.choose())

# rerun the chunk of the code until you get the m 
#m <- as.matrix(dtm)

#### change to upper case: e.g. required for some example: e.g. gene names

row.names(m)<- toupper(rownames(m))



##########
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)

set.seed(321)

wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))


## Ex:  make a word cloud of gene wiht somatic mutations over no. of patients##

#' 
#' 
#' ### 3. PCA plot - 2D
## ----fig.width = 6, fig.hei = 6------------------------------------------

library('corpcor')
library('ggplot2')

load('Data_Neve.RData')
Data <- dat.Neve
Info <- inf.Neve

dim(Data) 
# 19856 variables(genes), 50 samples 
dim(Info)
# head(Info) # cluster group

dx = as.matrix(Data)

# Method 1: Fast Singular Value Decomposition

svd = fast.svd(dx - rowMeans(dx)) # standardize by subject row Means

v = svd$v; rownames(v) = colnames(dx)


pca = list(v = v, d = svd$d)  

pcaVar = round((pca$d^2) / sum(pca$d^2) * 100, 2) # fraction of variace explained

x.lab = sprintf("PC1 (%.2f%% variance)", pcaVar[1])
y.lab = sprintf("PC2 (%.2f%% variance)", pcaVar[2])

# String Formatting, take the value of pcaVar to put in a string 

pcaDat = data.frame(sample=colnames(dx),
                      PC1=pca$v[,1], PC2=pca$v[,2],
                      group = Info$geneCluster)

pcaPlot = ggplot(pcaDat, aes(x=PC1, y=PC2, color=group)) +
    geom_point(stat="identity", size=3) +
    xlab(x.lab) + ylab(y.lab) +
    ggtitle(sprintf("PCA")) +
    theme(axis.ticks=element_blank(), axis.text.x=element_text(angle=-90))
pcaPlot

# show the 50 samples with PC1, PC2 are combination of 20K genes#
# remove color=group and try 

## PCA : Method  2##
# dim (t(dx)) #50*10986

GeneData<- data.frame( t(dx), grp=Info$geneCluster)
Gene.pca <- prcomp(GeneData[,1:19856], center = TRUE,scale. = F)
summary(Gene.pca) 


# Importance of components:
#   PC1      PC2      PC3      PC4      PC5
# Standard deviation     32.0867 23.08225 17.98058 16.65963 15.42048
# Proportion of Variance  0.1656  0.08572  0.05201  0.04465  0.03826
# Cumulative Proportion   0.1656  0.25136  0.30337  0.34802  0.38628

library(ggfortify)

autoplot(Gene.pca, data=GeneData, size=3)

# plot first PCA
autoplot(Gene.pca, data=GeneData, colour ='grp', size=3)
# label cell-lines
autoplot(Gene.pca, data=GeneData, colour ='grp',  shape = FALSE,label.size=3)

## Iris example
df <- iris[c(1, 2, 3, 4)]

summary(prcomp(df))
# Good amount of variance explained by first two PCAs

# unknown species class
autoplot(prcomp(df), data = iris)

# add color by species
autoplot(prcomp(df), data = iris, colour = 'Species')


#' 
#' 
#' ### 4. PCA plot - 3D
## ------------------------------------------------------------------------


#' 
#' - 4.1 Interactive plot  
library(rgl)

color <- factor(Info$geneCluster)
levels(color) <- c('cornflowerblue', 'coral3', 'green3')


plot3d(pca$v[,1:3],
       col = color, size = 10,
       xlab = 'PC1', ylab = 'PC2', zlab = 'PC3')
 
 # this opens a rgl window, can rotate
 # make it a bigger and can turn and rotate , then save
snapshot3d(file="persp3-new.png",fmt="png")



#' - 4.2  3D-plot  :static
#' 
## ----fig.width = 8, fig.height = 8---------------------------------------
library(scatterplot3d)

color <- factor(Info$geneCluster)
levels(color) <- c('cornflowerblue', 'coral3', 'green3')

angles <- c(60, 120, 240, 320) # preset some angles


scatterplot3d(x = pca$v[,1],y = pca$v[,2], z = pca$v[,3],
              main= paste("PCA 3D scatterplot,", "default angle"),
              xlab = 'PC1', ylab = 'PC2', zlab = 'PC3',
              color = color, pch = 16, cex.symbols = 2)
              

par(mfrow = c(2,2))

# with four user specified angle
for(i in 1: length(angles)){
  
  scatterplot3d(x = pca$v[,1],y = pca$v[,2], z = pca$v[,3],
              main= paste("PCA 3D,", "angle =", angles[i[]]),
              xlab = 'PC1', ylab = 'PC2', zlab = 'PC3',
              color = color, pch = 16, cex.symbols = 2,
              angle = angles[i])
 
}


#' 
#' ### 5. Venn diagram
## ------------------------------------------------------------------------

library(VennDiagram)
library(colorspace)

A <- sample(1:1000, 400, replace = FALSE)
B <- sample(1:1000, 600, replace = FALSE)
C <- sample(1:1000, 350, replace = FALSE)
D <- sample(1:1000, 550, replace = FALSE)

venn.diagram(list('A' = A, 'B' = B),
             filename = 'Venn-2.tiff',
             col = rainbow_hcl(2, alpha = 1),
             fill = rainbow_hcl(2, alpha = 0.5),
             cex = 3, cat.cex = 3)

venn.diagram(list('A' = A, 'B' = B, 'C' = C),
             filename = 'Venn-3.tiff',
             col = rainbow_hcl(3, alpha = 1),
             fill = rainbow_hcl(3, alpha = 0.5),
             cex = 3, cat.cex = 3)

venn.diagram(list('A' = A, 'B' = B, 'C' = C, 'D' = D),
             filename = 'Venn-4.tiff',
             col = rainbow_hcl(4, alpha = 1),
             fill = rainbow_hcl(4, alpha = 0.5),
             cex = 3, cat.cex = 3)

#' 
#' ### 6. Genomic plots - ggbio
#' 
## ------------------------------------------------------------------------
# source("https://bioconductor.org/biocLite.R"); biocLite("ggbio")

# BiocManager : install and manager packages from Bioconductor project #

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ggbio")
BiocManager::install("EnsDb.Hsapiens.v75")

library(ggplot2)
library(ggbio)
library(biovizBase)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75) #  annotation databases generated from Ensembl.


#### Circular plot: 22 chromosomes in one circle , also annotation


data("CRC", package = "biovizBase")


##.. single sample

gr.crc1 <- crc.gr[values(crc.gr)$individual == "CRC-1"]

p <- ggbio() +
 circle(gr.crc1, geom = "link", linked.to = "to.gr", aes(color = rearrangements)) +
 circle(gr.crc1, geom = "point", aes(y = score, size = tumreads), color = "red", grid = TRUE) + scale_size(range = c(1, 2.5)) +
 circle(mut.gr, geom = "rect", color = "steelblue") +
 circle(hg19sub, geom = "ideo", fill = "gray70") +
 circle(hg19sub, geom = "scale", size = 2) +
 circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
p


## break down the code to study the feature ##

#load CRC data from biovizBase
? CRC #a sample data of Genomic sequencing of colorectal adenocarcinomas

# loaded 3 data: [1] "crc.gr"  "hg19sub" "mut.gr"
# crc.gr
# GRanges object with 650 ranges and 17 metadata columns:

#   seqinfo: 22 sequences from an unspecified genome 

head(crc.gr,2 )
# contain sample, rarrangement, link info

names(values(crc.gr))
# [1] "individual"        "str1"              "class"             "span"             
# [5] "tumreads"          "normreads"         "gene1"             "gene2"            
# [9] "site1"             "site2"             "fusion"            "quality"          
# [13] "score"             "BPresult"          "validation_result" "to.gr"            
# [17] "rearrangements"   

head(hg19sub)
# seqnames      ranges strand
# <Rle>   <IRanges>  <Rle>
#   [1]        1 1-249250621      *
#   [2]        2 1-243199373      *
#   [3]        3 1-198022430      *
#   [4]        4 1-191154276      *
#   [5]        5 1-180915260      *
#   [6]        6 1-171115067      *
#GRanges object with 22 ranges and 0 metadata columns:
#seqinfo: 22 sequences from hg19 genome

head(mut.gr)
# GRanges object with 6 ranges and 10 metadata columns:


p1 <- ggbio()+
  circle(hg19sub, geom = "ideo", fill = "gray70"); p1 # Ideogram: graphic represent the chromosome

p1 <- p1+  circle(hg19sub, geom = "scale", size = 2) ; p1 #scale

p1 <- p1+   circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 4,col='red')
p1 #add labels

p1 <-  p1+circle(mut.gr, geom = "rect", color = "steelblue", radius=33); p1 
# add a "rectangle" track to show somatic mutation-looks like vertical segments.


#Next, add some "links" to show the rearrangement: intra and inter
#add a "point" track with grid background for rearrangement data and map ‘y‘ to variable
#      "score", map ‘size‘ to variable "tumreads", rescale the size to a proper size range.      

p1 <- p1 + circle(gr.crc1, geom = "point", aes(y = score, size = tumreads),
                  color = "red", grid = TRUE, radius = 30) +
                scale_size(range = c(1, 2.5)); p1  

#Finally, let’s add links and map color to rearrangement types. Remember you need to specify
#‘linked.to‘ parameter to the column that contain end point of the data.     

p1 <- p1 + circle(gr.crc1, geom = "link", linked.to = "to.gr", aes(color = rearrangements),
                  radius = 23)
p1
# try that in one step above , no need to set the radius ##



##.. multiple samples
# 9 single circular plots put together in one page, since we cannot keep too many
# tracks, we only keep ideogram and links. Here is one sample.

table(values(crc.gr)$individual)

#
#CRC-1 CRC-2 CRC-3 CRC-4 CRC-5 CRC-6 CRC-7 CRC-8 CRC-9 
#24     5   126    92    86    52    22    66   177 


# make a list to apply the same steps for each sample 
grl <- split(crc.gr, values(crc.gr)$individual)

## need "unit", load grid
library(grid)

crc.lst <- lapply(grl, function(gr.cur){
  print(unique(as.character(values(gr.cur)$individual)))
  cols <- RColorBrewer::brewer.pal(3, "Set2")[2:1]
  names(cols) <- c("interchromosomal", "intrachromosomal")
  p <- ggbio() + circle(gr.cur, geom = "link", linked.to = "to.gr",
  aes(color = rearrangements)) +
  circle(hg19sub, geom = "ideo",
  color = "gray70", fill = "gray70") +
  scale_color_manual(values = cols) +
  labs(title = (unique(values(gr.cur)$individual))) +
  theme(plot.margin = unit(rep(0, 4), "lines"))
 })

arrangeGrobByParsingLegend(crc.lst, widths = c(4, 1), legend.idx = 1, ncol = 3)

##############################################################################

#### Plot gene mode (ENSEMBL)

# In the example below we load an Ensembl based annotation package 
# for human gene and protein annotations defined in Ensembl version 75.
# It depends on  ensembldb package to retrieve gene/transcript/exons annotations 
# stored in an Ensembl based database 
# https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html

# Ensembl  transcript ID =  any spliced transcripts (ENST...) with overlapping coding sequence,
# 9 Ensembl transcript ID  shown here from GRCh37 genome

library(ensembldb )
ensdb <- EnsDb.Hsapiens.v75

(Tx <-transcripts(ensdb, filter = GeneNameFilter("GREB1")))

# Growth regulation by estrogen in breast cancer 1
# autoplot is a generic function to visualize various data object
#  which = A GRanges object to subset the result to plot the transcripts

autoplot(ensdb, which=GeneNameFilter("GREB1"))



#### Manhattan plot

# prep example data
snp <- read.table(system.file("extdata", "plink.assoc.sub.txt", package = "biovizBase"),
                  header = TRUE)
#head(snp,2)

# transform to GRanges object
gr.snp <- transformDfToGr(snp, seqnames = "CHR", start = "BP", width = 1)
gr.snp <- keepSeqlevels(gr.snp, as.character(1:22))

gr.snp <- gr.snp[!is.na(gr.snp$P)] 
values(gr.snp)$pvalue <- -log10(values(gr.snp)$P)

#head(gr.snp)

autoplot(gr.snp, geom = "point", coord = "genome", aes(y = pvalue))

plotGrandLinear(gr.snp, aes(y = pvalue), 
                color = c("#7fc97f", "#fdc086"),
cutoff = 3, cutoff.color = "blue", cutoff.size = 0.2)
# here cut 10e-3

# use different cut off/color
plotGrandLinear(gr.snp, aes(y = pvalue), color = c("#7fc97f", "#fdc086"),
                cutoff = 3.5, cutoff.color = "red", 
                cutoff.size = 0.2, linetype=2)

#####################################
#### survival plot 
#install.packages("survminer")
# check the cheat-sheet 

library(survminer)
library(survival)

# ? survival::lung 
# check censoring status
str(lung)
head(lung) #	censoring status 1=censored, 2=dead

# generate the estimate of survival distribution 
Lungfit <- survfit(Surv(time, status) ~ sex, data = lung)
plot(Lungfit, col=1:2)

ggsurvplot(Lungfit)

# customize plot area
ggsurvplot(Lungfit,  size = 1,  
           linetype = "strata", 
           risk.table = TRUE, fun='pct',
           risk.table.col = "strata", 
           break.time.by = 250, 
           xlab="Time in days",
           legend="bottom",
           legend.title="Sex",
           legend.labs=c("Male", "Female"),
           conf.int = TRUE, 
           pval = TRUE, 
           ggtheme=theme_gray() )

## check cumulative incidence (probability) of events ##


ggsurvplot(Lungfit,  size = 1,  
           linetype = "strata", 
           risk.table = TRUE, fun='event', ## change fun
           risk.table.col = "strata", 
           break.time.by = 100, 
           xlab="Time in days",
           legend="bottom",
           legend.title="Sex",
           legend.labs=c("Male", "Female"),
           conf.int = TRUE, 
           ggtheme=theme_bw() )

## make a forest plot of effects

lung$age <- ifelse(lung$age > 70, ">70","<= 70")
fit <- coxph( Surv(time, status) ~ sex + ph.ecog + age, data = lung)
summary(fit)

ggforest(fit, data=lung)











