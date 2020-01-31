##Usage: Rscript *.r vplot.txt fraction.txt bar-max-of-vplot y-max-of-geomline
args <- commandArgs(TRUE)
library("ggplot2")
library("reshape2")
library(Cairo)
library(RColorBrewer)
library(scales)
mytheme1 <- theme(plot.title = element_text(size=20, colour = "black",face = "bold"), axis.title = element_text(size=15,colour = "black",face = "bold"), axis.text = element_text(size=12,colour = "black",face = "bold"), axis.text.x = element_text(angle= 0,hjust=0.5, vjust=0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line = element_line(colour = "black", size = .5)) #,   panel.border = element_rect(colour = "black", fill = NA, size = 1)

mytheme2 <- theme_bw() + theme(plot.title = element_text(size=20, colour = "black",face = "bold"), axis.title = element_text(size=15,colour = "black",face = "bold"), axis.text = element_text(size=12,colour = "black",face = "bold"), axis.text.x = element_text(angle= 0,hjust=0.5, vjust=0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

m <- read.table(args[1], header=TRUE, sep="\t", as.is = TRUE, fill=TRUE, nrows=1000, comment.char = "",quote = "") #args[1]
#m <- m[,c(1,1752:2252)]
pk_size <- length(m[1,])-2

for (i in (-pk_size/2):(pk_size/2)) {names(m)[pk_size/2+2+i] <- i}
m1 <- melt(m, measure.vars = 2:(pk_size+2) , value.name = "CPB")
hp <- ggplot(m1, aes(x= as.numeric(as.character(variable)), y= fraglength, fill= CPB)) 
plotm <- hp + geom_tile() + scale_x_continuous(breaks = seq(-pk_size/2,pk_size/2,pk_size/10), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,200)) + scale_fill_gradientn(colours = c('white','red'), limits=c(0, as.numeric(args[3])), oob=squish) + labs(x = "Distance from Center (bp)", y = "Fragment Length (bp)") + mytheme1
mname <- paste(strsplit(args[1], 'txt', fixed=T)[[1]][1], "png", sep="")
ggsave(mname, plotm, width = 6.85, height = 4.97)

n <- read.table(args[2], header=TRUE, sep="\t", as.is = TRUE, fill=TRUE, nrows=1000, comment.char = "",quote = "")
#n <- n[c(175:225),]
names(n) <- c("distance", "80 bp", "100 bp", "120 bp", "140 bp", "168 bp")
nmelt <- melt(n, measure.vars = 2:6, value.name = "CPB")
lp <- ggplot(nmelt, aes(x= distance, y= CPB, colour= variable))
plotn <- lp + geom_line(size=1) + scale_x_continuous(expand = c(0, 0), breaks = seq(-pk_size/2,pk_size/2, pk_size/10)) + scale_y_continuous(expand = c(0, 0), limits = c(0, as.numeric(args[4]))) + labs(title = "",x = "Distance from Center (bp)", y = "CPB") + scale_colour_manual(name = " ", values = c("orangered", "orange", "limegreen", "deepskyblue", "magenta"), guide=F) + mytheme2 
nname <- paste(strsplit(args[2], 'txt', fixed=T)[[1]][1],"png", sep="")
ggsave(nname, plotn, width = 6.85, height = 4.97)




