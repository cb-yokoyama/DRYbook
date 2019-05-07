# written in 2013-11-22

### cpt.data()
### cpt.test()
### cpt.plot()


cpt.data <- function(i, seed = 123){
  set.seed(seed)
  data <- c (
    rnorm(i, mean = 0, sd = 1), # n = 100
    rnorm(i, mean = 1, sd = 1),
    rnorm(i, mean = 0, sd = 1),
    rnorm(i, mean = 0.2, sd = 1)
    )
  return(data)
}


cpt.test <- function(i, method = "PELT",penalty = "SIC", pen.value = 0,seed = 123) {
  require(changepoint)
  data <- cpt.data(i, seed = seed)
  # calculate changepoints
  m.cpt <- cpt.mean(data      = data,
                    method    = method,
                    penalty   = penalty,
                    pen.value = pen.value)
  return (m.cpt)
}

cpt.load <- function(data, method="PELT", penalty="SIC", pen.value=0){
  require(changepoint)
  data <- data
  m.cpt <- cpt.mean(data      = data,
                    method    = method,
                    penalty   = penalty,
                    pen.value = pen.value)
  return (m.cpt)  
}


cpt.plot <- function(cpt,
                     point.col = "grey",
                     cpt.col   = "#FF0000",
                     rect.col  = "#FF000040",
                     cpt.width = "4",
                     xlim      = NULL,
                     ylim      = c(0,1),
                     pos       = NULL,
                     axes     =FALSE,
                     main.title= NULL,ann=FALSE,cex=0.8,plot=TRUE){
  require(changepoint)
  # set display title param
  pen.value  <- round(pen.value(cpt), digits = 2)
  points     <- length(data.set(cpt))
  cpts       <- ncpts(cpt)
  nseg       <- cpts + 1
  # set title
  if(is.null(main.title)){
    main.title <- paste("pen.value = ", pen.value,
                        ", points = ", points,
                        ", nseg = ", nseg,
                        sep = "")
  }
  datay <- data.set.ts(cpt)
  if(length(pos) != 0){
    datax <- pos 
    xlab <- "pos"
  } else {
    datax <- index(data.set.ts(cpt))
    xlab <- "Index"
  }  
  # plotting
  if(plot==TRUE){
    plot(x = datax,
         y = datay,
         xlim = xlim,
         ylim = ylim,
         main = main.title,
         type = "p", 
         xlab = xlab,
         ylab = "% methylation",
         pch = 19, cex=cex, 
         axes=axes, ann=ann,
         col = point.col)
  } else {
    plot(x=1,y=1,
         xlim = xlim,
         ylim = ylim,
         main = main.title,
         type = "p", 
         xlab = xlab,
         ylab = "% methylation",
         pch = 19, cex=cex, 
         axes=axes, ann=ann,
         col = point.col)
  }
if(axes==TRUE){
  axis(1,labels=NA,lwd=2)
  axis(2,lwd=2)
}
  box(lwd=2)
  cpts <- c(0, cpt@cpts) 
  nseg <- ncpts(cpt)+1
  means <- param.est(cpt)$mean
  
    for(i in 1:nseg){
      segments(x0 = datax[cpts[i]+1], # cpts[1] => 0
               y0 = means[i],
               x1 = datax[cpts[i+1]],
               y1 = means[i],
               col= cpt.col, #ifelse(means[i]<50,cpt.col,"red")
               lwd=cpt.width)
        rect(xleft=datax[cpts[i]+1],
           ybottom=-10,
            xright=datax[cpts[i+1]],
             ytop=110,
             col=rect.col,
             lty="dotted",border=NA)
    }
}   

cpt.ans <- function(i){
segments(x0=seq(from=1,to=3*i+1,by=i),
        y0=c(0,1,0,0.2),
        x1=seq(from=i,to=(3+1)*i,by=i),
        y1=c(0,1,0,0.2),
         col="red",lwd=4,lty="dashed"
        )  
}

grid.cpt.plot <- function(cpt,
                     point.col = "grey",
                     cpt.col   = "#FF0000",
                     rect.col  = "#FF000040",
                     cpt.width = "4",
                     xlim      = NULL,
                     pos       = NULL,
                     main.title= NULL,...){
  require(changepoint)
  # set display title param
  pen.value  <- round(pen.value(cpt), digits = 2)
  points     <- length(data.set(cpt))
  cpts       <- ncpts(cpt)
  nseg       <- cpts + 1
  # set title
  if(is.null(main.title)){
    main.title <- paste("pen.value = ", pen.value,
                        ", points = ", points,
                        ", nseg = ", nseg,
                        sep = "")
  }
  datay <- data.set.ts(cpt)
  if(length(pos) != 0){
    datax <- pos 
    xlab <- "pos"
  } else {
    datax <- index(data.set.ts(cpt))
    xlab <- "Index"
  }  
  # plotting
  grid.points(x = datax,
       y = datay,
  #     xlim = xlim,
  #     main = main.title,
  #     type = "p", 
  #     xlab = xlab,
  #     ylab = "% methylation",
       pch = 19,
  #     axes=T,
       gp=gpar(col = point.col),...)

  cpts <- c(0, cpt@cpts) 
  nseg <- ncpts(cpt)+1
  means <- param.est(cpt)$mean
  
    for(i in 1:nseg){
      grid.segments(x0 = datax[cpts[i]+1], # cpts[1] => 0
               y0 = means[i],
               x1 = datax[cpts[i+1]],
               y1 = means[i],
               gp=gpar(col= cpt.col, #ifelse(means[i]<50,cpt.col,"red")
               lwd=cpt.width))
    }
}   

# filtering methylation context & threshold
filterContCov <- function(data,class="CG",threshold=5){
  
  # context
  data_cg <- subset(data, {
      class==class
      }, select=c(assembly,position,mc,h) )
  
  # threshold  
  data_m5p <- subset(data, {
    h>=threshold &  class==class
    }, select=c(assembly,position,mc,h) )

  
  if(is.unsorted(data$position)){
    data_m5ps <- data_m5p[order(data_m5p$position),]
  } else {
    data_m5ps <- data_m5p
  }
  
  data_m5pw <- within(data_m5ps,{
    seqid <- assembly
    pos <- position  
    mratio <- round(mc/h,digits=2)
    rm(mc,h,position)
  })[,c("seqid","pos","mratio")]
  return(data_m5pw)
}


cpt.plot2 <- function(cpt,
                     point.col = "grey",
                     cpt.col   = "#FF0000",
                     rect.col  = "#FF000040",
                     cpt.width = "4",
                     xlim      = NULL,
                     ylim      = c(0,1),
                     cpg       = NULL,
                     axes     =FALSE,
                     main.title= NULL,ann=FALSE,cex=0.8){
  #require(changepoint)
  # set display title param
  #pen.value  <- round(pen.value(cpt), digits = 2)
  #points     <- length(data.set(cpt))
  #cpts       <- ncpts(cpt)
  cpts <- nrow(cpt) - 1
  nseg       <- cpts + 1 
  # set title
  if(is.null(main.title)){
#    main.title <- paste("pen.value = ", pen.value,
 #                       ", points = ", points,
  #                      ", nseg = ", nseg,
   #                     sep = "")
  }
  #datay <- data.set.ts(cpt) func2
  datay <- cpg[,"mratio"]
  if(nrow(cpg) != 0){
    datax <- cpg[,"pos"] 
    xlab <- "pos"
  } else {
    datax <- index(data.set.ts(cpt))
    xlab <- "Index"
  }  
  # plotting
  plot(x = datax,
       y = datay,
       xlim = xlim,
       ylim = ylim,
       main = main.title,
       type = "p", 
       xlab = xlab,
       ylab = "% methylation",
       pch = 19, cex=cex, 
       axes=axes, ann=ann,
       col = point.col)
  if(axes==TRUE){
    axis(1,labels=NA,lwd=2)
    axis(2,lwd=2)
  }
  box(lwd=2)
  #cpts <- c(0, cpt[,"end"]) 
  #nseg <- ncpts(cpt)+1
  #means <- param.est(cpt)$mean
  
  for(i in 1:nseg){
    segments(x0 = cpt[i,"start"], # cpts[1] => 0
             y0 = cpt[i,"mratio"],
             x1 = cpt[i,"end"],
             y1 = cpt[i,"mratio"],
             col= cpt.col, #ifelse(means[i]<50,cpt.col,"red")
             lwd=cpt.width)
  }
}   
