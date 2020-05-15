## Draw a heart in R ##

dat<- data.frame(t=seq(0, 2*pi, by=0.1) )
xhrt <- function(t) 16*sin(t)^3
yhrt <- function(t) 13*cos(t)-5*cos(2*t)-2*cos(3*t)-cos(4*t)

dat$y=yhrt(dat$t)
dat$x=xhrt(dat$t)

with(dat, plot(x,y, type="l", axes=F, xlab="", ylab=""))
with(dat, polygon(x,y, col="pink"))

text(mean(dat$x ), mean(dat$y),  "Happy Valentine's Day",col='red',cex=2.5, font=2)  
text(mean( dat$x ), mean(dat$y)-3,  Sys.Date() ,col='red',cex=1.5, font=3) 
 # add names and make a card 