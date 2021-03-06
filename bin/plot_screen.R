####################################
#Read in file for gentrap analysis #
####################################

read.screen <- function(filename, header=TRUE) {
  gentrap.Pval <- read.table(filename, header=header, sep="\t")
  gentrap.Pval <- gentrap.Pval[,1:7]
  screen.header <- c("genesymbol",
                     "num_insertions_screen",
                     "no_insertions_screen",
                     "num_insertions_control",
                     "no_insertions_control",
                     "pval",
                     "adj_pval")
  colnames(gentrap.Pval) = screen.header
  gentrap.Pval$adj_pval <- as.numeric(as.character(gentrap.Pval$adj_pval))
  return(gentrap.Pval)
}

###########################
#Read in Refseq file      #
###########################
#TODO: reformat refseq file, remove empty lines with "strange" characters

read.all.refseq <- function(inputfile2, header=TRUE){
  descTab <- read.delim(file=inputfile2, header=header, 
                        sep="\t", dec=".", stringsAsFactors=FALSE)
  descTab <- descTab[,1:ncol(descTab)-1]  # remove last, empty column
  if (!header){
    names(descTab) <- c("refseq_id","chrom","strand","start","end",
                        "CDS_start","CDS_end","exon_count","genesymbol", 
                        "gene_description")
  }
  return(descTab)
}

###########################################
#assign chromosomal position to gene names#
###########################################

assign.chr.pos <- function(gentrap.screens.mult, refseq.genes){
  # merge tables
  resultTab <- merge(gentrap.screens.mult, refseq.genes, 
                     by="genesymbol", all.x=TRUE, sort=FALSE)
  return(resultTab)
}

##################
#load and process#
##################
load.process.data <-function(filename, inputRefgene){
  
    gentrap.screens <- read.screen(filename, header=TRUE)
    refseq.genes <- read.all.refseq(inputRefgene, header=FALSE)
    wo_descTab <- subset(refseq.genes, start!="") #remove "strange" characters
    wo_descTab$start <- as.numeric(wo_descTab$start)
    #merge files on gene_symbol
    gentrap.screens.mult <- assign.chr.pos(gentrap.screens, wo_descTab)
    return(gentrap.screens.mult)   
}

subset.genes <-function(screen, num.genes=NULL, siglevel=NULL){
  gene.list <- screen
  #subset gene list 
  if(!is.null(num.genes)){
    gene.list <- gene.list[order(gene.list$adj_pval),]$genesymbol[1:num.genes] 
    gene.list <- as.character(gene.list)
    
  }
  else {
    gene.list <- gene.list[gene.list$adj_pval<siglevel,]$genesymbol
    gene.list <- as.character(gene.list)
  }
  
  return(gene.list)
}

#############
#Creat plot#
#############
create.bubble.plot <- function(genetrap.screens, selection.list, upper.limit=0.1, 
                               lower.limit=0.01, breaks, limits, max_s, 
                               scale_text_size=TRUE){
  
  #convert pvalues to -log10(pvalue)
  genetrap.screens$mlog10Pval <- -log10(genetrap.screens$adj_pval)
  
  # label x-axis
  plot.col = c("Chromosomal position of genes")
  x.axis.lab = c("start")
  
  # little hack, needed to substitute upper.limit[1] with 
  # value in geom_hline
  genetrap.screens$upper.limit <- upper.limit
  genetrap.screens$lower.limit <- lower.limit
  
  #correct x-axis label for chromosomal position
  ticks=NULL
  lastbase=0
  genetrap.screens$pos=0 
    
  #order by chromosome and position
  chrOrder <- c(paste("chr",1:22,sep=""),"chrX") #chrY not included why??
  genetrap.screens$chrom <- factor(genetrap.screens$chrom, levels=chrOrder)
  genetrap.screens <- genetrap.screens[order(genetrap.screens$chrom, 
                                             genetrap.screens$start),]
  genetrap.screens <- na.omit(genetrap.screens) #remove NA values 
    
  for (i in 1:length(unique(genetrap.screens$chrom))) {
    chrom.list <- unique(genetrap.screens$chrom)
    if (chrom.list[i]=="chr1") {
      genetrap.screens[genetrap.screens$chrom=="chr1", ]$pos = genetrap.screens[genetrap.screens$chrom=="chr1", ]$start  
    }
      
    else {
      lastbase=lastbase+tail(subset(genetrap.screens,
                                    genetrap.screens$chrom==chrom.list[i-1])$start, 1)
      genetrap.screens[genetrap.screens$chrom==chrom.list[i], ]$pos <- genetrap.screens[genetrap.screens$chrom==chrom.list[i], ]$start+lastbase
    }
    ticks=c(ticks, genetrap.screens[genetrap.screens$chrom==chrom.list[i], ]$pos[floor(length(genetrap.screens[genetrap.screens$chrom==chrom.list[i], ]$pos))])
  }   
  
  #colorblind friendly colors
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                 "#999999", "#E69F00", "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                 "#999999", "#E69F00", "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", "#D55E00")
  
  #creat ggplot
  g_obj <- ggplot(data=genetrap.screens, aes(x=pos, y=mlog10Pval, size=num_insertions_screen))
    g_obj = g_obj + geom_point(aes(fill=chrom), colour="black", pch=21, alpha=1/4)
    g_obj = g_obj + geom_point(data=genetrap.screens[(genetrap.screens$genesymbol %in% selection.list),], 
                               aes(x=pos, y=mlog10Pval, fill=chrom), colour="black", pch=21)
    g_obj = g_obj + scale_fill_manual(values=cbPalette, name = "chromosome \ncolor")
    g_obj = g_obj + scale_size_area(max_size = max_s, name = "number of \ninsertions")
    if (scale_text_size) g_obj = g_obj + geom_text(data = genetrap.screens[(genetrap.screens$genesymbol %in% selection.list),], 
                                                   aes(x=pos, y=mlog10Pval, label=genesymbol), colour="black", vjust=-0.5, hjust=1.1)
    else g_obj = g_obj + geom_text(data = genetrap.screens[(genetrap.screens$genesymbol %in% selection.list),], 
                                   aes(x=pos, y=mlog10Pval, label=genesymbol), size=10, colour="black", vjust=-0.5, hjust=1.1)
    g_obj = g_obj + geom_text(data = genetrap.screens[(genetrap.screens$genesymbol %in% selection.list),], 
                              aes(x=pos, y=mlog10Pval, label=num_insertions_screen), size=4, colour="black", vjust=1.5, hjust=-0.6)
    g_obj = g_obj + scale_x_continuous(name=plot.col[1], breaks=ticks, labels=(unique(genetrap.screens$chrom)))
    g_obj = g_obj + scale_y_log10(name=expression(paste("Significance [-", log[10], "(",italic("FDR-corrected p-value"),")]")),
                                  breaks=breaks, limits=limits)
    g_obj = g_obj + theme_bw()
    g_obj = g_obj + theme(axis.text.x=element_text(angle = -90, hjust = 0, size=15), 
                          axis.text.y=element_text(colour="black",size=15), 
                          panel.grid.major=element_line(colour = "grey55", size = 0.2), 
                          axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))
    g_obj = g_obj + guides(fill = guide_legend(ncol=2))
    g_obj = g_obj + geom_hline(aes(yintercept=-log10(upper.limit[1])), colour="#56B4E9", linetype="dashed")
    g_obj = g_obj + geom_hline(aes(yintercept=-log10(lower.limit[1])), colour="#990000", linetype="dashed")
    
  return(g_obj)
} 

optimize_limits <- function(adj_pval){
  
  y_high <- -log10(min(adj_pval)) # lowest pval
  y_low <- -log10(max(adj_pval)) # highest pval
  
  breaks = NULL
  break_p = NULL
  seq_a = (-10:10)
  
  for(i in seq_a){
    break_p <- 10^(i)
    if (break_p > y_low){
      seq_b <- (i-1):10
      for(i in seq_b){
        break_p <- 10^(i)
        breaks <- cbind(breaks, break_p)
        if (break_p > y_high){break}
      }
    break
    }
  }
    
  return(as.vector(breaks))
}
  
print_plot <- function(filename, g_obj){
  pdf(file=out_f, width=14, height=10, useDingbats=FALSE)
  print(g_obj)
  dev.off()
}


library("ggplot2")
options <- commandArgs(trailingOnly=TRUE)
refseq_f = options[1]
in_f = options[2]
out_f = options[3]

#refseq_f = "~/Dropbox/src_code/bioinf_workflows/data/All_genes_data_ct_srt_sed.txt"
#in_f = "~/Dropbox/src_code/bioinf_workflows/data/TEST.fisher_test.txt"
#out_f = "~/Dropbox/src_code/bioinf_workflows/data/TEST.bubble_plot.pdf"

screen <- load.process.data(in_f, 
                            refseq_f)
signif_genes <- subset.genes(screen, 
                             num.genes=NULL, 
                             siglevel=0.05)
breaksforplot <- optimize_limits(screen$adj_pval)
print(breaksforplot)
g_obj <- create.bubble.plot(screen, signif_genes, 
                            upper.limit=0.05,
                            lower.limit=0.05,
                            breaks=breaksforplot,
                            limits=c(breaksforplot[1], 
                                     breaksforplot[length(breaksforplot)]),
                            max_s=30,
                            scale_text_size=FALSE)
print_plot(out_f, 
           g_obj)
