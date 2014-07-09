Bubble plots
============

```R
create.bubble.plot <- function(genetrap.screens, selection.list, 
                               upper.limit=0.1, lower.limit=0.01, 
                               breaks, limits, range, scale_text_size=TRUE, 
                               stat="BH"){
  
    #convert pvalues to -log10(pvalue)
    genetrap.screens$mlog10Pval <- -log10(genetrap.screens[[stat]]);
  
    #label x-axis
    plot.col = c("Chromosomal position of genes")
    x.axis.lab = c("start")
    
    #little hack, needed to substitute upper.limit[1] with value in geom_hline
    genetrap.screens$upper.limit <- upper.limit
    genetrap.screens$lower.limit <- lower.limit
  
    #correct x-axis label for chromosomal position
    ticks=NULL
    lastbase=0
    genetrap.screens$pos=0 
    
    #order by chromosome and position
    chrOrder <- c(paste("chr",1:22,sep=""),"chrX") #chrY not included why??
    genetrap.screens$chrom <- factor(genetrap.screens$chrom, levels=chrOrder)
    genetrap.screens <- genetrap.screens[order(genetrap.screens$chrom, genetrap.screens$start),]
    genetrap.screens <- na.omit(genetrap.screens) #remove NA values 
    
    for (i in 1:length(unique(genetrap.screens$chrom))) {
      chrom.list <- unique(genetrap.screens$chrom)
      if (chrom.list[i]=="chr1") {
        genetrap.screens[genetrap.screens$chrom=="chr1", ]$pos = genetrap.screens[genetrap.screens$chrom=="chr1", ]$start  
      }
      
      else {
        lastbase=lastbase+tail(subset(genetrap.screens,genetrap.screens$chrom==chrom.list[i-1])$start, 1)
        genetrap.screens[genetrap.screens$chrom==chrom.list[i], ]$pos <- genetrap.screens[genetrap.screens$chrom==chrom.list[i], ]$start+lastbase
      }
      ticks=c(ticks, genetrap.screens[genetrap.screens$chrom==chrom.list[i], ]$pos[floor(length(genetrap.screens[genetrap.screens$chrom==chrom.list[i], ]$pos))])
    }   
  
  #colorblind friendly colors
  cbPalette <- c("#999999", "#E69F00", 
                 "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", 
                 "#D55E00", "#CC79A7",
                 "#999999", "#E69F00", 
                 "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", 
                 "#D55E00", "#CC79A7",
                 "#999999", "#E69F00", 
                 "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", 
                 "#D55E00")
  
  #creat ggplot
  g_obj <- ggplot(data=genetrap.screens, aes(x=pos, y=mlog10Pval, size=insertions))
    g_obj = g_obj + geom_point(aes(fill=chrom), colour="black", pch=21, alpha=1/4)
    g_obj = g_obj + geom_point(data=genetrap.screens[(genetrap.screens$gene %in% selection.list),], aes(x=pos, y=mlog10Pval, fill=chrom), colour="black", pch=21)
    g_obj = g_obj + scale_fill_manual(values=cbPalette, name = "chromosome \ncolor")
    g_obj = g_obj + scale_size(range = range, name = "number of \ninsertions")
    if (scale_text_size) g_obj = g_obj + geom_text(data = genetrap.screens[(genetrap.screens$gene %in% selection.list),], aes(x=pos, y=mlog10Pval, label=gene), colour="black", vjust=-0.5, hjust=1.1)
    else g_obj = g_obj + geom_text(data = genetrap.screens[(genetrap.screens$gene %in% selection.list),], aes(x=pos, y=mlog10Pval, label=gene), size=10, colour="black", vjust=-0.5, hjust=1.1)
    g_obj = g_obj + geom_text(data = genetrap.screens[(genetrap.screens$gene %in% selection.list),], aes(x=pos, y=mlog10Pval, label=insertions), size=4, colour="black", vjust=1.5, hjust=-0.6)
    g_obj = g_obj + scale_x_continuous(name=plot.col[1], breaks=ticks, labels=(unique(genetrap.screens$chrom)))
    g_obj = g_obj + scale_y_log10(name=expression(paste("Significance [-", log[10], "(",italic("FDR-corrected p-value"),")]")), breaks=breaks, limits=limits)
    g_obj = g_obj + theme_bw()
    g_obj = g_obj + theme(axis.text.x=element_text(angle = -90, hjust = 0, size=15), axis.text.y=element_text(colour="black",size=15), panel.grid.major=element_line(colour = "grey55", size = 0.2), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))
    g_obj = g_obj + guides(fill = guide_legend(ncol=2))
    g_obj = g_obj + geom_hline(aes(yintercept=-log10(upper.limit[1])), colour="#56B4E9", linetype="dashed")
    g_obj = g_obj + geom_hline(aes(yintercept=-log10(lower.limit[1])), colour="#990000", linetype="dashed")
    
    return(g_obj)
}
