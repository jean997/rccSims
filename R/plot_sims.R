#'Plot coverage of example simulations
#'@param R A list produced by \code{example_sim}
#'@param simnames Names of CIs to plot. The order they are given will be the plotting order.
#'@param legend.names Names to print in legend
#'@param cols,shapes colors and shapes to use.
#'@param ltys Line types.
#'@param main Title
#'@param proportion Proportion of statistics to output
#'@param legend.position Where to put the legend.
#'@param y.axis.off Don't label the y axis?
#'@param span. Value of span to pass to loess. If present, coverage will be smoothed.
#'@return A ggplot object
#'@export
plot_coverage <- function(R, simnames, cols, shapes, ltys, 
                          legend.names=NULL, main="", proportion=0.2,
                          y.axis.off=FALSE, y.range=c(0, 1),
                          legend.position=c(0.28, 0.4), span=NULL){
  which.keep <- which(R$simnames %in% simnames)
  if(is.null(legend.names)) legend.names <- simnames
  ncis <- length(which.keep)
  p <- dim(R$COVERAGE)[2]
  k <- ceiling(proportion*p)
  avg.coverage <- matrix(nrow=k, ncol=ncis)
  nms <- c()
  for(i in 1:ncis){
    w <- which.keep[i]
    nms <- c(nms, R$simnames[w])
    avg.coverage[,i] <- rowMeans(R$COVERAGE[w, , ])[1:k]
  }
  avg.coverage <- data.frame(avg.coverage)
  names(avg.coverage) <- nms
  avg.coverage$Rank <- 1:k
  avg.coverage.long <- gather(avg.coverage, "Method", "RCC", -Rank)
  if(is.null(y.range)) y.range =range(avg.coverage.long$RCC)
  #Re-order factor levels
  avg.coverage.long$Method <- factor( as.character(avg.coverage.long$Method),
                                      levels=simnames)
  m <- avg.coverage.long$Method
  o <- c()
  for(i in 1:length(simnames)){
    o <- c(o, which(m==simnames[i]))
  }
  avg.coverage.long <- avg.coverage.long[o, ]
  if(is.null(span)){
    h <- ggplot(avg.coverage.long, aes(x=Rank)) + geom_hline(aes(yintercept=0.9)) +
      geom_point(aes(y=RCC,  color=Method, shape=Method)) +
      scale_shape_manual(values=shapes, labels=legend.names)
  }else{
    h <- ggplot(avg.coverage.long, aes(x=Rank)) + geom_hline(aes(yintercept=0.9)) +
      geom_line(aes(y=RCC,  color=Method, linetype=Method),
                span=span, stat="smooth", method="loess", level=0, alpha=0.8, size=1.5) +
      scale_linetype_manual(values=ltys, labels=legend.names)
  }
  h <- h+  theme_bw(base_size = 14) + ylim(y.range) +
    scale_color_manual(values=cols, labels=legend.names) +
    theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())

  h <- h + theme(legend.position=legend.position,
                      legend.key = element_rect(color="white", fill="white", size = 0.6, linetype=0),
                      legend.background = element_rect(fill="transparent"),
                      legend.title= element_text(size=15),
                      legend.text=element_text(size=9))
  if(y.axis.off) h <- h + theme(axis.title.y = element_blank())
  if(!main == "") h <- h + ggtitle(main) + theme(plot.title=element_text(size=20))
  return(h)

}

#'Plot average width of example simulations
#'@param R A list produced by \code{example_sim}
#'@param simnames Names of CIs to plot. The order they are given will be the plotting order.
#'@param legend.names Names to print in legend
#'@param cols,shapes colors and shapes to use.
#'@param main Title
#'@param proportion Proportion of statistics to output
#'@param legend.position Where to put the legend.
#'@param y.axis.off Don't label the y axis?
#'@param span. Value of span to pass to loess. If present, coverage will be smoothed.
#'@return A ggplot object
#'@export
plot_width <- function(R, simnames, cols, shapes, ltys, legend.names=NULL,
                               main="", proportion=0.2, span=NULL,
                               y.axis.off=FALSE, y.max=NULL,
                               legend.position="none"){
  which.keep <- which(R$simnames %in% simnames)
  if(is.null(legend.names)) legend.names <- simnames
  p <- dim(R$WIDTH)[2]
  k <- ceiling(proportion*p)
  ncis <- length(which.keep)
  ylim=range(apply(R$WIDTH[which.keep, , ], MARGIN=1,  FUN=function(x){ z <- rowMeans(x); return(z[1:k])}), na.rm=TRUE)
  if(!is.null(y.max)) ylim[2] <- y.max
  avg.width <- matrix(nrow=k, ncol=ncis)
  nms <- c()
  for(i in 1:ncis){
    w <- which.keep[i]
    nms <- c(nms, R$simnames[w])
    if(!is.null(span)){
      y <- rowMeans(R$WIDTH[w, , ])[1:k]
      x <- 1:k
      f <- loess(y~x, span=span)
      avg.width[f$x,i] <- f$fitted
    }else{
      avg.width[,i] <- rowMeans(R$WIDTH[w, , ])[1:k]
    }
  }
  avg.width <- data.frame(avg.width)
  names(avg.width) <- nms
  avg.width$Rank <- 1:k
  
  avg.width.long <- gather(avg.width, "Method", "Average Width", -Rank)
  avg.width.long$Method <- factor( as.character(avg.width.long$Method),
                                   levels=simnames)
  m <- avg.width.long$Method
  o <- c()
  for(i in 1:length(simnames)){
    o <- c(o, which(m==simnames[i]))
  }
  avg.width.long <- avg.width.long[o, ]

  if(is.null(span)){
    h <- ggplot(avg.width.long, aes(x=Rank)) +
      geom_point(aes(y=`Average Width`, group=Method,  color=Method, shape=Method)) +
      scale_shape_manual(values=shapes, labels=legend.names)
  }else{
    h <- ggplot(avg.width.long, aes(x=Rank)) +
      geom_line(aes(y=`Average Width`, group=Method,  color=Method, linetype=Method),
               alpha=0.8, size=1.5)+
      scale_linetype_manual(values=ltys, labels=legend.names)
  }
  h <- h +  theme_bw(base_size = 14) +
    scale_color_manual(values=cols, labels=legend.names) +
    theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())

  h <- h + theme(legend.position=legend.position,
                 legend.key = element_rect(color="white", fill="white", size = 0.6, linetype=0),
                 legend.background = element_rect(fill="transparent"),
                 legend.title= element_text(size=15),
                 legend.text=element_text(size=9))

  if(y.axis.off) h <- h + theme(axis.title.y = element_blank())
  if(!main == "") h <- h + ggtitle(main) + theme(plot.title=element_text(size=20))
  h <- h + ylim(ylim)
  return(h)

}


make_sim_legend <- function(legend.names, cols, ltys){
  n <- length(legend.names)
  #points <- data.frame(x=rep(1, 5), y = rev(seq(1, 3, length.out=5)))
  points <- data.frame(y=rep(1, n), x = seq(1, 2.25, length.out=n))
  
  dist <- (2.25 -1)/n
  points$left = points$x-dist*0.4
  points$right = points$x + dist*0.4
  
  points$labs <- legend.names
  points$lty <-ltys
  points$cols <- cols
  
  p <-ggplot(points) +
    geom_segment(aes(x=left, xend=right, y=y, yend=y), lty=points$lty,
                 colour=points$cols, lwd=1.3)+
        annotate(geom="text", y=rep(1.1, n), x=points$x,
             label=points$labs, size=4) +
    xlim(0.8, 2.5) + ylim(0.97, 1.2) +
    theme_bw() + theme(panel.grid=element_blank(), axis.title=element_blank(),
                       axis.text=element_blank(),
                       panel.border=element_blank(), axis.ticks=element_blank())
  
  #ggsave(p, file="~/Dropbox/cfdr-jean/for_AOAS/img/sim_legend.png", height=1, width=8, units="in", dpi=300)
  return(list("plot"=p, "info"=points))
}


plot_cis <- function(rank, ci, truth, prop=1, plot.truth=FALSE){
  dat <- data.frame("Rank"=rank, "ciL"=ci[,1], "ciU"=ci[,2], "truth"=truth)
  dat$covered <-dat$ciL <= truth & dat$ciU >= truth
  if(prop < 1){
    nkeep <- ceiling(nrow(dat)*prop)
    dat <- dat[dat$Rank <=nkeep,]
  }
  p <- ggplot(dat) + geom_linerange(aes(x=Rank, ymin=ciL, ymax=ciU, col=covered)) 
  if(plot.truth) p <- p + geom_point(aes(x=Rank, y=truth, col=covered))
  p <- p + 
    scale_color_discrete(name="Parameter\nCovered") + 
    theme_bw()+
    theme(panel.grid=element_blank())
  return(p)
}


