# Library import
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
all_genes <- genes(txdb, filter=c(exon_chrom="chr1"))

# DMR data import(liftOvered)
dmr <- read.csv("hg19_DMR_output.bed",sep="\t",head=F)

names(dmr)<- c("chromosome","start","end","id")
dmr.1 <- dmr[which(dmr$chromosome=="chr1"),c("start","end")]
dmr.1$str <- "*"
dmr.1$color <- topo.colors(16)[4]

source("DataImport.R")
source("R/cpttest.R")
source("R/plotOperons.R")
start <- dmr.1$start-15000

#pngfileNames <- paste0("image/Fig2B_chr1_",start,"_",start+30000,".png",sep="")

#for(i in 1:length(start)){
i <- 6 # for FOXD2
xCoord <- c(start[i]-2500,start[i]+35000-2500)
  
Trackwidth <- 0.125
oldpar <- par()

# epsfileNames <- paste0("image/Fig2B_chr1_",xCoord[1],"_",xCoord[2],".eps",sep="")
# pngfileNames <- paste0("image/Fig2B_chr1_",xCoord[1],"_",xCoord[2],".png",sep="")
# 
# plot.new()
# png(pngfileNames[i],width = 5.2, height=6.4,units = "in",pointsize=11,res=72)
#postscript(epsfileNames,colormodel="cmyk",width=5.2,height=6.4,pointsize=11,family="ArialMT")
  my.cex=0.7

  Pen.ceiling <- 1.0- 0.15

CGI<- as.data.frame(rbind(
  c(47899126,47899398),
  c(47909713,47911020),
  c(47902794,47905518),
  c(47881897,47883065),
  c(47915640,47915952),
  c(47899662,47900385)
))
names(CGI)<-c("start","end")

CGI$str<- "*"
CGI$color<-"green"

# AdditonalTracks
# - CpG island
# - Gene Annotation
# - penalty panel
query <- GRanges(seqnames="chr1",ranges=IRanges(xCoord[1],xCoord[2]),strand="*")
ann <- as.data.frame(all_genes)[subjectHits(findOverlaps(query,all_genes)),]

Ann.ceiling <- 1.0

# Track 0a; Annnotation Track
par(fig=c(0,1,Ann.ceiling-0.075,Ann.ceiling-0), 
    mar = c(bottom=0.2, left=5, top=0.2, right=0.2),xaxs="i")
if(nrow(ann)!=0){
  
  ann$color <- "gray"
  # ann$color<-ifelse(ann$strand=="+","blue","red")
  ann$str <- "*" #ifelse(ann$strand=="+",1,-1)
  plotOperons(ann[,c("start","end","str","color", "gene_id")], xCoord[1], xCoord[2], 
              axes=F,xlab="",lwd=2)
  box(lwd=2)
  mtext("Gene", side=2, line=0.5, cex.lab=1,las=1, col="black")
}

# Track 0b: CpG Island
par(fig=c(0,1,Ann.ceiling-0.15,Ann.ceiling - 0.075), 
    mar = c(bottom=0.2, left=5, top=0.2, right=0.2),xaxs="i",new=T)
plotOperons(CGI[,c("start","end","str","color")], xCoord[1], xCoord[2], heightOfOperon = 1.5,
            axes=F,xlab="",lwd=2)
box(lwd=2)
mtext("CpG island", side=2, line=0.5, cex.lab=1,las=1, col="black")


# Track 1: pv18
par(fig=c(0,1,Pen.ceiling-Trackwidth*1,Pen.ceiling),
    mar = c(bottom=0.2, left=5, top=0.2, right=0.2),xaxs="i",new=T)
cpt.plot(cpt.H1_chr1.pv18,pos=H1_chr1[,"pos"],main.title="",axes=F, point.col=gray(0.3),
         rect.col="#00FF0000",xlim=xCoord,ann=FALSE,cex=my.cex)
mtext("1.8", side=2, line=0.5, cex.lab=1,las=1, col="black")
axis(2,at=seq(0.2,0.8,0.2),label=NA,tck=1,lty="dashed",col=gray(0.5))

# Track 2: pv14
par(fig=c(0,1,Pen.ceiling-Trackwidth*2,Pen.ceiling-Trackwidth*1), mar = c(bottom=0.2, left=5, top=0.2, right=0.2),xaxs="i",new=T)
cpt.plot(cpt.H1_chr1.pv14,pos=H1_chr1[,"pos"],main.title="",axes=F, point.col=gray(0.3),
         rect.col="#00FF0000",xlim=xCoord,cex=my.cex)
mtext("1.4", side=2, line=0.5, cex.lab=1,las=1, col="black")
axis(2,at=seq(0.2,0.8,0.2),label=NA,tck=1,lty="dashed",col=gray(0.5))

# Track 3: pv10
par(fig = c(0,1,Pen.ceiling-Trackwidth*3,Pen.ceiling-Trackwidth*2), 
    mar = c(bottom=0.2, left=5, top=0.2, right=0.2), xaxs="i",new=T)
cpt.plot(cpt.H1_chr1.pv10,pos=H1_chr1[,"pos"],main.title="",axes=F, point.col=gray(0.3),
         rect.col="#00FF0000",xlim=xCoord,cex=my.cex)
mtext("1.0", side=2, line=0.5, cex.lab=1,las=1, col="black")
axis(2,at=seq(0.2,0.8,0.2),label=NA,tck=1,lty="dashed",col=gray(0.5))

# Track 4: pv06
par(fig = c(0,1,Pen.ceiling-Trackwidth*4,Pen.ceiling-Trackwidth*3), 
    mar = c(bottom=0.2, left=5, top=0.2, right=0.2),xaxs="i",new=T)
cpt.plot(cpt.H1_chr1.pv06,pos=H1_chr1[,"pos"],main.title="",axes=F, point.col=gray(0.3),
         rect.col="#00FF0000",xlim=xCoord,cex=my.cex)
mtext("0.6", side=2, line=0.5, cex.lab=1,las=1, col="black")
axis(2,at=seq(0.2,0.8,0.2),label=NA,tck=1,lty="dashed",col=gray(0.5))

# Track 5: pv02
par(fig = c(0,1,Pen.ceiling-Trackwidth*5,Pen.ceiling-Trackwidth*4),
    mar = c(bottom=0.2, left=5, top=0.2, right=0.2),xaxs="i",new=T)
cpt.plot(cpt.H1_chr1.pv02,pos=H1_chr1[,"pos"],main.title="",axes=F, point.col=gray(0.3),
         rect.col="#00FF0000",xlim=xCoord,cex=my.cex)
mtext("0.2", side=2, line=0.5, cex.lab=1,las=1, col="black")
axis(2,at=seq(0.2,0.8,0.2),label=NA,tck=1,lty="dashed",col=gray(0.5))
   
par(fig=c(0,1,0.1,Pen.ceiling-Trackwidth*5),
    mar = c(bottom=0.2, left=5, top=0.2, right=0.2),xaxs="i",new=T)

cpt.plot(cpt.H1_chr1.pv12,pos=H1_chr1[,"pos"],main.title="",axes=F, point.col="#00FF0000",cpt.col="#00FF0000",
         rect.col="#00FF0000",xlim=xCoord,cex=my.cex,ylim=c(0,2.4),plot=F)
my.lwd="3"
xpos_pv02=H1_chr1[cpt.H1_chr1.pv02@cpts,"pos"]
xpos_pv04=H1_chr1[cpt.H1_chr1.pv04@cpts,"pos"]
xpos_pv06=H1_chr1[cpt.H1_chr1.pv06@cpts,"pos"]
xpos_pv08=H1_chr1[cpt.H1_chr1.pv08@cpts,"pos"]
xpos_pv10=H1_chr1[cpt.H1_chr1.pv10@cpts,"pos"]
xpos_pv12=H1_chr1[cpt.H1_chr1.pv12@cpts,"pos"]
xpos_pv14=H1_chr1[cpt.H1_chr1.pv14@cpts,"pos"]
xpos_pv16=H1_chr1[cpt.H1_chr1.pv16@cpts,"pos"]
xpos_pv18=H1_chr1[cpt.H1_chr1.pv18@cpts,"pos"]
xpos_pv20=H1_chr1[cpt.H1_chr1.pv20@cpts,"pos"]

xpos <- xpos_pv02[!(xpos_pv02 %in% xpos_pv04)]
points(x = xpos, y=rep(0.2,length(xpos)),col="red",pch=19)

xpos <- xpos_pv04[!(xpos_pv04 %in% xpos_pv06)]
points(x = xpos, y=rep(0.4,length(xpos)),col="red",pch=19)

xpos <- xpos_pv06[!(xpos_pv06 %in% xpos_pv08)]
points(x = xpos, y=rep(0.6,length(xpos)),col="red",pch=19)

xpos <- xpos_pv08[!(xpos_pv08 %in% xpos_pv10)]
points(x = xpos, y=rep(0.8,length(xpos)),col="red",pch=19)

xpos <- xpos_pv10[!(xpos_pv10 %in% xpos_pv12)]
points(x = xpos, y=rep(1.0,length(xpos)),col="red",pch=19)

xpos <- xpos_pv12[!(xpos_pv12 %in% xpos_pv14)]
points(x = xpos, y=rep(1.2,length(xpos)),col="red",pch=19)

xpos <- xpos_pv14[!(xpos_pv14 %in% xpos_pv16)]
points(x = xpos, y=rep(1.4,length(xpos)),col="red",pch=19)

xpos <- xpos_pv16[!(xpos_pv16 %in% xpos_pv18)]
points(x = xpos, y=rep(1.6,length(xpos)),col="red",pch=19)

xpos <- xpos_pv18[!(xpos_pv18 %in% xpos_pv20)]
points(x = xpos, y=rep(1.8,length(xpos)),col="red",pch=19)

xpos <- xpos_pv20
points(x = xpos, y=rep(2.0,length(xpos)),col="red",pch=19)

abline(h=seq(0.2,1.8,0.4),lty="dashed")
axis(2, at=c(0,1.0,2.0),labels=c(0,1.0,2.0),las=1)
axis(1)

par(oldpar)

