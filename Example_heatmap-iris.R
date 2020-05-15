#' ---
#' title: "Using color in R"
#' output: word_document
#' ---
#' 
## ----setup, cache=FALSE, include=FALSE, warning = FALSE------------------

library(ggplot2)
library(grid)
library(gridExtra)

# R color Palettes packages
library(RColorBrewer)
library(colorspace)
library(colorRamps)

# clustering for color branches
library(dendextend)
# heatmap 
library(gplots)
library(hexbin)
library(MASS)


#' 
#' - 3.4 Exercises: 
#' 
#' - 3.4.1 Correlation heatmaps with gplots (use R data iris):  
#' 1. Use colorRampPalette to generate color palette (list colors and how many shades of colors): colorRampPalette(collect.colors)(number.shades)
#' 2. Heatmap 1: map heatmap colors to positive/negative correlation only;
#' 3. Heatmap 2: map heatmap colors to ranges of correlation coefficient;
#' 4. Heatmap 3: add cell borders (use: colsep/rowsep argument) and cell text (use cellnote argument)
#' 5. Heatmap 4: set missing data for diagonal (use: diag(matrix), omit dendrograms and legend key)
#' 
#' 
## ----fig.width=10, fig.height=10-----------------------------------------
data(iris)
head(iris)
cormat <- cor(iris[,1:4], method='spearman', use='pairwise.complete.obs')
cormat2<- round(cormat,2)
is.matrix(cormat2)
#' 
# heatmap 1: 4 choose 2= 6 numbers, sysmmetric matrix

#cormat2<- as.matrix(cormat2)

#Heatmap1
heatmap.2(cormat, trace='none',  
          col = colorpanel(2, low='darkblue',  high = 'brown3') , 
          margins =c(10,10), main="Heatmap 1")

#Heatmap2
heatmap.2(cormat2, trace='none',  
          col = colorpanel(50, low='darkblue', mid = 'yellow',
             high = 'darkred'),  margins =c(12,12), sepcolor="black",
          main="Heatmap 2")

#Heatmap3
heatmap.2(cormat2, trace='none',  
          col = colorpanel(200, low='darkgreen', mid = 'yellow',
                                                  high = 'red'),  
          cellnote=cormat2, notecol="black", notecex=1.2, margins =c(12,12),
       sepwidth=c(0.01,0.01), rowsep = 0:4, colsep=0:4, main="Heatmap 3",
        sepcolor="black")

#Heatmap4
cormat3 <- cormat2
diag(cormat3)<-NA

heatmap.2(cormat3, trace='none', col = colorpanel(200, low='darkgreen', mid = 'yellow',
                                                   high = 'red'),  
          cellnote=cormat2, dendrogram='none',
          notecol="black", notecex=1.2, margins =c(12,12), 
          key = F,main="Heatmap 4" )









