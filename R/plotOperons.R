library(org.Hs.eg.db)

convertToOperon <- function(x1, x2, strand, heightOfOperon=1, ratioOfTriangle=0.3, y1=1){
  # if strand =="1":
  #
  # (x1,y2)      (xAngle,y2)  
  #     ________________
  #    |                \
  # _ _| _ _ _ _ _ _ _ _ \ _ _ _
  #    |                 / (x2,y3)
  #    |________________/  
  # (x1,y1)            
  #
  y2 <- y1 + heightOfOperon
  y3 <- y1 + heightOfOperon/2
  if(strand == 1){
    xAngle <- (1-ratioOfTriangle)*(x2-x1)+x1
    x <- c(x1, x1, xAngle, x2, xAngle)
    y <- c(y1, y2, y2,     y3, y1)
    return(list(x, y))
  }else if(strand == -1){
    xAngle <- ratioOfTriangle*(x2-x1)+x1
    x <- c(x1, xAngle, x2, x2, xAngle)
    y <- c(y3, y2, y2, y1, y1)
    return(list(x, y))
  } else {
    x <- c(x1, x1, x2, x2)
    y <- c(y1, y2, y2, y1)
    return(list(x, y))
  }
}

plotOperons <- function(dat, start, end, axes=TRUE, line=FALSE,
                        heightOfOperon=0.7, 
                        ratioOfTriangle=0.3, y1=1, lwd=2, xlab="[kb]"){
  genome <- c(start, end)
  height <- c(0,3)
  # Define plot space
  plot(genome, height, type = "n", xlab = xlab, ylab="", axes=F)
  if(axes==TRUE){
    axis(1)
  }
  # Draw horizontal line
  if(line==TRUE){
    lines(c(start, end), c(y1+heightOfOperon/2, y1+heightOfOperon/2), lwd=lwd)
  }
  if(is.factor(dat[,4])){
    dat[,4] <- as.character(dat[,4])
  }
  for(i in 1:nrow(dat)){
    coo <- convertToOperon(dat[i,1], dat[i,2], dat[i,3], 
                           heightOfOperon=heightOfOperon, 
                           ratioOfTriangle=ratioOfTriangle,
                           y1=y1)
    polygon(coo[[1]], coo[[2]], col=dat[i,4], lwd=lwd)
    if(dat$str != "*"){
      offset <- ifelse(dat$str == "+", dat$start[i], dat$start[i])
      text(offset, 2.5, get(dat$gene_id[i], org.Hs.egSYMBOL), adj=0, cex=0.7)
    }
  }
}

