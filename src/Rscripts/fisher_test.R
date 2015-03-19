## Merge count tables by genesymbol and perform fisher-test &
## correct for multiple testing, FDR

read.input.control <- function(input){
  file1 <- read.table(file=input, sep="\t")
  colnames(file1) <- c("genesymbol", "num_insertions")
  total_insertions <- sum(file1$num_insertions)
  file1$no_insertions <- total_insertions - file1$num_insertions
  return(file1)
}

read.input.screen <- function(input){
  file1 <- read.table(file=input, sep="\t", header=TRUE)
  return(file1)
}

merge.tables <- function(screen, control){
  # merge files by genesymbol
  merged_table <- merge(screen, control, all=TRUE, by="genesymbol")
  # only use complete cases
  complete_cases <- merged_table[complete.cases(merged_table),]
  return(complete_cases)
}

cal.fisher <- function(merged_table){
  
.cal.fisher.test <- function(df){
  mat <- matrix(as.numeric(df[c(2:5)]), ncol=2)
  f <- fisher.test(as.table(mat), alternative="greater")
  return(c(df[1], df[2], df[3], df[4], df[5], f$p.value))
}
  row.names(merged_table) <- merged_table$genesymbol
  fisher_results <- apply(merged_table, 1, .cal.fisher.test)
  trans.fisher.results <- t(fisher_results)
  trans.fisher.results <- as.data.frame(trans.fisher.results, stringsAsFactors=FALSE)
  colnames(trans.fisher.results) <- c("genesymbol", "num_insertions.x", 
                                      "no_insertions.x", "num_insertions.y", 
                                      "no_insertions.y", "pval")
  trans.fisher.results$pval <- as.numeric(trans.fisher.results$pval)
  sorted <- trans.fisher.results[order(trans.fisher.results$pval),]
  # correct for multiple testing, benjamini-hochberg, fdr-correction
  sorted$adj <- p.adjust(as.numeric(sorted[,6]), "fdr")
  return(sorted)  
}

write2file <- function(fisher_table, filename){
  fn_out <- filename
  write.table(fisher_table, file=paste(fn_out),
              append=FALSE, quote=FALSE, 
              row.names=FALSE, col.names=TRUE, 
              dec=".", sep="\t")
}

options <- commandArgs(trailingOnly=TRUE)
#screen_f = options[1]
#control_f = options[2]
#file_out = options[3]
screen_f = "~/Dropbox/src_code/bioinf_workflows/data/TEST.count_table.txt"
control_f = "~/Dropbox/src_code/bioinf_workflows/data/control_set_gene+screen_incorrect_count.txt"
file_out = "~/Dropbox/src_code/bioinf_workflows/test/TEST.fisher_test.txt"

screen <- read.input.screen(screen_f)
control <- read.input.control(control_f)
merged_t <- merge.tables(screen, control)
fisher_t <- cal.fisher(merged_t)
write2file(fisher_t, file_out)
