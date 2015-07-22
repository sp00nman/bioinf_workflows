library(stringr)




parse_gtf <- function(inputDir, inputFile, outDir)
{  # load file 
   print(paste("Loading of", inputFile, ".."))
   setwd(inputDir)
   annotTab0 <- read.delim(file=inputFile, header=FALSE, sep="\t", dec=".", stringsAsFactors=FALSE)
   names(annotTab0) <- c("chrom","source","feature","start","end","score","strand","frame","attributes")
   
   # get ids, save in separate columns
   print(paste("Parse ids .."))
   idSource <- str_split(annotTab0$attributes, " ")
   vIds <- str_replace(unlist(lapply(idSource, "[[", 2)), ";", "")
   vSymbols <- str_replace(unlist(lapply(idSource, "[[", 4)), ";", "")
   vTc <- str_replace(unlist(lapply(idSource, "[[", 6)), ";", "")
   annotTab1 <- transform(annotTab0, gene_id=vIds, gene_symbol=vSymbols, tc_id=vTc)
   
   # get exons, start and stop codons
   annotTab2 <- annotTab1[grep("on", annotTab1$feature),]
   
   # create final table with selected columns
   finalTab <- annotTab2[,c(1,3,4,5,7,10:12)]
   names(finalTab) <- c("chrom", "feature", "start", "end", "strand", "gene_id", "gene_symbol", "transcript_id")
   print(head(finalTab))
   print(nrow(finalTab))
   
   # save table
   outfile <- paste0(str_split(inputFile, "\\.")[[1]][1], "_parsed.txt")
   setwd(outDir)
   print(paste("Saving table .."))
   write.table(finalTab, file=outfile, append=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)
   print(paste(outfile, "saved in", outDir))  
 
}



args <- commandArgs(trailingOnly=TRUE)
print(args)
if(length(args)<3)
{  print("Rscripts parse_gtf_for_gene_ins_plots_150717.R INPUT_DIR GTF_FILE OUT_DIR")
   stop()
}
inputDir <- args[1]
inputFile <- args[2]
outDir <- args[3]

inputDir <- "Z:/groups/lab_kralovics/fschischlik/projects/prj_diploid_screen/annotation/"
inputFile <- "ucsc_hg19_ensembl_73_genes.gtf"
outDir <- "C:/Data/KBM7/test_files"
parse_gtf(inputDir, inputFile, outDir)



