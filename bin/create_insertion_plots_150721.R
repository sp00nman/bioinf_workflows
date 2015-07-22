# Author: Doris Chen

# REQUIRED PACKAGES
library(grid)

# FUNCTIONS
create_insertion_plots <- function(annotFilePath, hitFilePath, insFilePath, outDir, fdrCutoff, minDistFactor, screenName, plotOption, plotWidth, plotHeight)
{  # defaults
   padjIndex <- 7  # in hits table
   
   # create output folders if not existing
   if(!file.exists(outDir))
   {  dir.create(outDir)
      print(paste(outDir,"created"))
   }
     
   # for plots
   fileNameBase <- paste0(screenName,"_mdf",minDistFactor)
   color_blue <- rgb(0,0,200,150,maxColorValue=255)
   color_red <- rgb(200,0,0,150,maxColorValue=255)
   triangleHeight <- 2
      
   # get screen result 
   print(paste("Loading of gene hit table", hitFilePath, "..."))
   hitTab <- read.delim(file=hitFilePath, header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)
   print(head(hitTab))
   print(paste(nrow(hitTab), "hits found."))
   
   # get gene list
   selListIndex <- which(hitTab[,padjIndex] < fdrCutoff)
   geneList <- hitTab[selListIndex,geneIndex]
   print(geneList)   
   print(paste(length(geneList)," genes with adjusted p-value <",fdrCutoff)) 
      
   # load insertion position table
   print(paste("Loading of insertions table",insFilePath,"..."))
   insTab <- read.delim(file=insFilePath, header=FALSE, sep="\t", dec=".", stringsAsFactors=FALSE)
   names(insTab) <- c("chrom", "start", "end", "id", "length", "strand", "cigar")
   print(head(insTab))
   print(paste(nrow(insTab), "insertions found."))
         
   # load annotation files
   print(paste("Loading of annotation table",annotFilePath,"..."))
   featureTab <- read.delim(file=annotFilePath, header=TRUE, sep="\t", dec=".", stringsAsFactors=FALSE)
   print(head(featureTab))
   print(nrow(featureTab))
   
   # go through gene list, create plots with insertions (one per isoform)
   for(currGene in geneList)
   {  print("")
      print("")
      print(currGene)
      
      # get transcripts
      tcList <- unique(featureTab$transcript_id[which(featureTab$gene_symbol==currGene)])
      print(paste(length(tcList),"transcripts found."))
        
      for(currTc in tcList)
      {  print("")
         print(paste0(currGene,", ",currTc))
         
         # get features from subTab
         subTabFeatures <- featureTab[which(featureTab$transcript_id==currTc),]
         subTabExon <- subTabFeatures[grep("exon", subTabFeatures$feature),]
         currChrom <- unique(subTabExon$chrom)
         currStrand <- unique(subTabExon$strand)
         vStart <- subTabExon$start
         vEnd <- subTabExon$end
         xRange <- c(min(vStart), max(vEnd))
         currStart <- subTabFeatures[which(subTabFeatures$feature=="start_codon"), ifelse(currStrand=="-", 4, 3)][1]   # in case of multiple start codons, first one picked
         currStop <-  subTabFeatures[which(subTabFeatures$feature=="stop_codon"), ifelse(currStrand=="-", 3, 4)][1]  # in case of multiple stop codons, first one picked
         if(is.na(currStart))  # in case of NR
         {  currStart <- min(subTabFeatures$start)
            currStop <- max(subTabFeatures$end)
         } 
         
         # get insertions from subTab
         insSubTab <- insTab[which(insTab$chrom==currChrom & insTab$start>=xRange[1] & insTab$start<=xRange[2]),]
         insSubTab <- insSubTab[order(insSubTab$start),]
         vInsSense <- insSubTab[which(insSubTab$strand==currStrand),]
         print(paste(nrow(vInsSense),"sense insertions found."))
         vInsAntisense <- insSubTab[which(insSubTab$strand!=currStrand),]
         print(paste(nrow(vInsAntisense),"antisense insertions found."))
         
         if(nrow(insSubTab)>0)
         {  # plot
            setwd(outDir)
            if(plotOption=="png")
            {  plotFileName <- paste0(fileNameBase, "_", currGene, "-", currTc, ".png")  
               png(file=plotFileName, width=plotWidth, height=plotHeight, res=300) 
            } else
            if(plotOption=="pdf")
            {  plotFileName <- paste0(fileNameBase, "_", currGene, "-", currTc, ".pdf")   
               pdf(file=plotFileName, width=plotWidth, height=plotHeight)
            } 
            
            # preparation of plotting region
            grid.newpage()
            pushViewport(plotViewport(c(5, 4, 2, 2)))
            pushViewport(dataViewport(xRange, seq(-80*triangleHeight*0.7,20), name="plotRegion"))  
            grid.xaxis(gp=gpar(fontsize=16))
            
            # chrom. name
            grid.text(currChrom, x=unit(vStart[1], "native"), y=0.03, just=c("left","top"), gp=gpar(col="black", fontsize=17))
            
            # gene name
            grid.text(paste0(currGene, " (",currTc,")"), x=unit(vStart[1], "native"), y=0.99, just=c("left","top"), gp=gpar(col="black", fontsize=25))
            
            # introns
            vAll <- c(vStart[-1],vEnd[-(length(vEnd))])
            vX <- vAll[order(vAll)]
            vY <- rep(0, length(vX))   
            for(i in seq(1, length(vAll), 2))
            {  grid.lines(x=unit(vX[i:(i+1)], "native"), y=unit(vY[i:(i+1)], "native"), name="intron", gp=gpar(col="grey", lty="solid", lwd=3))
            }
            
            # orientation arrows   
            if(length(currStop)>0 & length(currStart)>0)
            {  arrowHeight <- triangleHeight/2
               arrowWidth <- round((xRange[2]-xRange[1])/250)
               arrowStep <- 3000*sign(currStop-currStart)
               if(currStrand=="+")
               {  for(arrStart in seq(currStart+20,currStop-20,arrowStep))
                  {  grid.lines(x=c(arrStart, arrStart+arrowWidth, arrStart), y=c(-arrowHeight, 0, arrowHeight), default.units="native", gp=gpar(col="grey", lty="solid", lwd=2), name="orientation")
                  }
               }
               if(currStrand=="-")
               {  for(arrStart in seq(currStart-20,currStop+20,arrowStep))
                  {  grid.lines(x=c(arrStart, arrStart-arrowWidth, arrStart), y=c(-arrowHeight, 0, arrowHeight), default.units="native", gp=gpar(col="grey", lty="solid", lwd=2), name="orientation")
                  }
               }
            }
            
            # exons
            vX <- vStart
            vY <- rep(0, length(vX))
            vWidth <- vEnd - vX
            exonHeight <- 5   
            # x and y specify left, bottom point (according to just)
            grid.rect(x=unit(vX, "native"), y=unit(vY, "native"), width=unit(vWidth,"native"), height=unit(exonHeight,"native"), name="exon", just=c("left","center"), gp=gpar(col="black", lty="solid", fill="black"))
            
            # insertions
            triangleWidth <- round((xRange[2]-xRange[1])/150)
            # sense
            vDistance <- vInsSense$start[-1] - vInsSense$start[-(length(vInsSense$start))]   
            insYOri <- -4
            minY <- insYOri
            insY <- insYOri
            for(i in seq_along(vInsSense$start))
            {  if(currStrand=="+")
               {  pointer <- vInsSense$start[i]+triangleWidth
               } else   
               if(currStrand=="-")
               {  pointer <- vInsSense$start[i]-triangleWidth
               }
               grid.polygon(x=c(vInsSense$start[i], vInsSense$start[i], pointer, vInsSense$start[i]), y=c(insY-triangleHeight/2, insY+triangleHeight/2, insY, insY-triangleHeight/2), default.units="native", gp=gpar(col=color_red, fill=color_red, lty="solid", lwd=1), name="insertion_sense")
               if( vDistance[i]<triangleWidth/minDistFactor & i<length(vInsSense$start))
               {  insY <- insY - triangleHeight
                  minY <- min(insY, minY)
               }
               else
                insY <- insYOri
            }
            # antisense
            vDistance <- vInsAntisense$start[-1] - vInsAntisense$start[-(length(vInsAntisense$start))]
            #insYOri <- minY - triangleHeight 
            insY <- insYOri 
            for(i in seq_along(vInsAntisense$start))
            {  if(currStrand=="+")
               {  pointer <- vInsAntisense$start[i]-triangleWidth
               } else   
               if(currStrand=="-")
               {  pointer <- vInsAntisense$start[i]+triangleWidth 
               }
               grid.polygon(x=c(vInsAntisense$start[i], vInsAntisense$start[i], pointer, vInsAntisense$start[i]), y=c(insY-triangleHeight/2, insY+triangleHeight/2, insY, insY-triangleHeight/2), default.units="native", gp=gpar(col=color_blue, fill=color_blue, lty="solid", lwd=1), name="insertion_sense")
               if( vDistance[i]<triangleWidth/minDistFactor & i<length(vInsAntisense$start))
               {  insY <- insY - triangleHeight
               } else
                 insY <- insYOri
            }         
            
            if(plotOption=="png" | plotOption=="pdf")
              dev.off()
         } else
         {  print("No plot created.")          
         }
      }
   }
}


# USER INPUT
args <- commandArgs(trailingOnly=TRUE)
print(args)
if(length(args)<7)
{  print("Rscripts create_insertion_plots_150721.R ANNOTATION_FILE(gtf format) HITS_FILE INSERTIONS_FILE OUTPUT_DIR(will be created if not existing) FDR_CUTOFF(e.g. 0.05) SCREEN_NAME(for plot file names) PLOT_FORMAT(png|pdf) PLOT_WIDTH(e.g. 3000|24) PLOT_HEIGHT(e.g. 2400|16) [MIN_DIST_FACTOR(integer; default=10)]")
   stop()
}
annotFilePath <- args[1]
hitFilePath <- args[2]
insFilePath <- args[3]  # containing screen results
outDir <- args[4]  # for processed table
fdrCutoff <- as.numeric(args[5])
screenName <- args[6]
plotOption <- args[7]
plotWidth <- as.numeric(args[8])
plotHeight <- as.numeric(args[9])

if(length(args)==10)
{  minDistFactor <- as.integer(args[10])  # any number except 0; if 10 -> if distance between triangles is smaller than 1/10th of the triangle width, the second triangle is plotted next to (and not below) the first
} else
{  minDistFactor <- 10  # any number except 0; if 10 -> if distance between triangles is smaller than 1/10th of the triangle width, the second triangle is plotted next to (and not below) the first (thus the higher the mdf, the more densely the insertions are plotted)
}

# annotFilePath <- "C:/Data/KBM7/test_files/ucsc_hg19_ensembl_73_genes_parsed.txt"
# hitFilePath <- "C:/Data/KBM7/test_files/TEST_bsf.fisher_test.txt"
# insFilePath <- "C:/Data/KBM7/test_files/TEST_bsf.fix_header_rm2bp_insertions.bed"  # containing screen results
# outDir <- "C:/Data/KBM7/test_files/gene_insertion_plots"  # for processed table
# fdrCutoff <- 0.001
# screenName <- "TEST"
# plotOption <- "png"
# plotWidth <- 3000
# plotHeight <- 2400
# minDistFactor <- 10
create_insertion_plots(annotFilePath, hitFilePath, insFilePath, outDir, fdrCutoff, minDistFactor, screenName, plotOption, plotWidth, plotHeight)




