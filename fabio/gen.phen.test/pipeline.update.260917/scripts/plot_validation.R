#!/usr/bin/env Rscript
# Standalone scientific summaries; base R only. Usage:
# Rscript scripts/plot_validation.R summaries/RUN/matched/strategy_comparison.tsv figures [functional_dir]
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<2) stop('Supply strategy_comparison.tsv and output directory')
x <- read.delim(args[1], check.names=FALSE, stringsAsFactors=FALSE)
x$smoke <- tolower(as.character(x$smoke)) == 'true'
dir.create(args[2], recursive=TRUE, showWarnings=FALSE)
order.ids <- c('P0','N0','N1','N2','N3','N4','N5','R0','R1','P1','N6','P2')
x <- x[order(match(x$strategy_id,order.ids), x$hypothesis_id),]
save.figure <- function(name, draw, width=9, height=5) {
  pdf(file.path(args[2],paste0(name,'.pdf')), width=width,height=height,useDingbats=FALSE)
  draw(); dev.off()
  png(file.path(args[2],paste0(name,'.png')),width=width,height=height,units='in',res=300)
  draw(); dev.off()
}
plot.discoveries <- function() {
  par(mar=c(5.5,5,2,1),las=1,cex=1.1)
  if (!any(is.finite(x$common_query_fraction))) {
    plot.new(); text(.5,.5,'No shared completed, testable genes\nFractions are not estimable')
    return(invisible(NULL))
  }
  groups <- order.ids[order.ids %in% x$strategy_id]
  values <- lapply(groups,function(g) x$common_query_fraction[x$strategy_id==g])
  names(values) <- groups
  boxplot(values,ylim=c(0,max(.05,x$common_query_fraction,na.rm=TRUE)*1.12),ylab='Query / common testable background',xlab='',
          col=ifelse(groups=='P0','#983548',ifelse(groups %in% c('R0','R1'),'#c7d4dc','#91a9b8')),outline=FALSE)
  for(i in seq_along(values)) points(rep(i,length(values[[i]])),values[[i]],pch=16,cex=.7,col='#34424c66')
  mtext('Each point is one hypothesis; non-completed genes excluded',side=1,line=3.3,cex=.85)
  if(any(x$smoke)) mtext('SYNTHETIC SOFTWARE TEST — NOT BIOLOGICAL RESULTS',side=3,col='#983548',cex=.8)
}
save.figure('01_query_background_fraction',plot.discoveries)
fixed <- x[!x$strategy_id %in% c('R0','R1'),]
if(nrow(fixed)) save.figure('02_positions_before_after_filter',function(){
  par(mar=c(5.5,5,2,1),las=1,cex=1.1)
  bx <- barplot(fixed$significant_positions,names.arg=fixed$strategy_id,col='#c2d2dc',border=NA,
                ylab='Nominally significant primary positions',ylim=c(0,max(1,fixed$significant_positions)*1.2))
  rect(bx-.5,0,bx+.5,fixed$retained_positions,col='#435d70',border=NA)
  legend('topright',c('Before dense-cluster pruning','Retained'),fill=c('#c2d2dc','#435d70'),bty='n',cex=.8)
  mtext('P1 / N6 use 15 linked cycles: compare internally, not as 100-cycle arms',side=1,line=3.3,cex=.75)
  if(any(x$smoke)) mtext('SYNTHETIC SOFTWARE TEST',side=3,col='#983548',cex=.8)
})
if(length(args)>=3 && file.exists(file.path(args[3],'conditional_tail_fractions.tsv'))) {
  tails <- read.delim(file.path(args[3],'conditional_tail_fractions.tsv'),stringsAsFactors=FALSE)
  if(nrow(tails)) save.figure('03_conditional_benchmark',function(){
    par(mar=c(5,5,2,1),las=1)
    plot(seq_len(nrow(tails)),tails$conditional_tail_fraction,log='y',xaxt='n',pch=19,col='#983548',
         xlab='',ylab='Conditional upper-tail fraction',ylim=range(c(pmax(1e-5,tails$monte_carlo_exceedance_ci_lower),tails$monte_carlo_exceedance_ci_upper)))
    axis(1,seq_len(nrow(tails)),tails$null_family)
    arrows(seq_len(nrow(tails)),pmax(1e-5,tails$monte_carlo_exceedance_ci_lower),seq_len(nrow(tails)),tails$monte_carlo_exceedance_ci_upper,angle=90,code=3,length=.08)
    mtext('Intervals describe Monte Carlo exceedance uncertainty; not a universal association test',side=1,line=3,cex=.7)
    if(any(x$smoke)) mtext('SYNTHETIC SOFTWARE TEST',side=3,col='#983548',cex=.8)
  })
}
