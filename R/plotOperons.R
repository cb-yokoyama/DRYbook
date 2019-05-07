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
    #  print(dat[i, 4])
      polygon(coo[[1]], coo[[2]], col=dat[i,4], lwd=lwd)
    }
  }

# dat <- data.frame(start=c(5, 50, 70),
#                   end=c(30, 60, 90),
#                   strand=c(1, -1, 1),
#                   color=c("orange", "green", "magenta"))
# 
# layout(matrix(c(0,1,0,2,3,3), 3, 2, byrow = TRUE), 
#        widths=c(1,3), heights=c(2,2,1))
# layout.show(3)
# oldpar<-par()
# par(mfrow=c(3,1),mar = c(0.2, 0.2, 0.2, 0.2))
# cpt.plot(cpt.aic.h1,pos=mc_h1_1f[,"pos"],main.title="",
#          rect.col="#00FF0000",xlim=xCoord)
# cpt.plot(cpt.aic.h1,pos=mc_h1_1f[,"pos"],main.title="",
#          rect.col="#00FF0000",xlim=xCoord)
# plotOperons(dat, 0, 100)
# par(oldpar())
