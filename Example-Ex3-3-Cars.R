#' ---
#' title: "R Graphics Devices"
#' output: word_document
#' ---
#' 
## ----setup, cache=FALSE, include=FALSE, warning = FALSE------------------
data(mtcars)
head(mtcars)
str(mtcars)
?mtcars

### File1 ##
? pdf()

pdf('Cars_plot1.pdf', width = 10, height = 8) # default 7 inches
 # plot1 : hp:	 Gross horsepower; ranked high to low with labels

 n<- nrow(mtcars)
 mtcars2<- mtcars[rev(order(mtcars$hp)),]

 plot(1:n, sort(mtcars$hp), type='n' , ylim=c(0, 350), xlab="rank", ylab="gross horsepower", main="Data Obtained in mid-1970s")
  for (i in 1:n)
  {arrows(i,0, i,mtcars2$hp[i], length=0, lwd=2 , col=1)
   text(i+0.3,mtcars2$hp[i]-25, srt=90, paste( row.names(mtcars2)[i]), cex=0.8, col=4 )}

 
 #plot2
  par(mar=c(5,8,4,3))
  mtcars2<- mtcars[order(mtcars$hp),]
  counts<- as.matrix(t(mtcars2$cyl+  mtcars2$carb))
 
  A<-barplot(counts, horiz=TRUE, col='gray90', xlab="count", ylab="", main="Data ordered by gross horsepower" )
  axis(2, at=A,  label= row.names(mtcars2), col.axis=1, cex.axis=0.9, las=2, tick = F, mgp=c(3,0.02,0), pos=-0.2)
  barplot(as.matrix(t(mtcars2$cyl)), horiz=TRUE, col='gray20', xlab="", ylab="", main="", add=T )
  legend('bottomright',   fill=c('gray20', 'gray90' ), legend=c('cylinders', 'carburetors') )
  
dev.off()
  
### File 2##

ggplot(mtcars, aes(x=carb, y=qsec)) + 
  geom_jitter(size=2)+
  stat_smooth(method="lm", col='red', se=FALSE)+
  labs(title='Cars data in mid-1970s', x="Carb", y="Qsec")

ggsave("Cars_plot2.tiff", w=4, h=4, dpi = 300 , units='in')

