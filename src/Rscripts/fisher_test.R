## Fisher's exact test

get_fisher <- function(df){
  mat <- matrix(as.numeric(df[c(1:4)]), ncol=2)
  f <- fisher.test(as.table(mat), alternative="greater")
  return(c(df[1], df[2], df[3], df[4], f$p.value))
}

path="path"
fn1="f1"
file1 <- read.table(file=paste(path,fn1,sep=""), sep="\t")
colnames(file1) <- c("number_insertions","gene_symbol")

fn2="f2"
file2 <- read.table(file=paste(path,fn2,sep=""), sep="\t")
colnames(file2) <- c("number_insertions","gene_symbol")

# merge files by gene symbol
# leave non-common gene names as 0
merged.table <- merge(file1, file2, all=TRUE, by="gene_symbol")
merged.table[is.na(merged.table)] <- 0
total.x <- colSums(merged.table[2:3])[1]-merged.table$number.insertions.x
total.y <- colSums(merged.table[2:3])[2]-merged.table$number.insertions.y

merged.table.counts <- cbind(number.insertions.y=merged.table$number.insertions.y,
                             total.y,
                             number.insertions.x=merged.table$number.insertions.x, 
                             total.x)
rownames(merged.table.counts) <- merged.table$gene_symbol

fisher_results <- apply(merged.table.counts, 1,  get_fisher)
trans.fisher.results <- t(fisher_results)
colnames(trans.fisher.results) <- c("number_insertions.y", "total.y", "number_insertions.x", "total.x","pval")
trans.fisher.results <- as.data.frame(trans.fisher.results)
sorted <- trans.fisher.results[order(trans.fisher.results$pval),]

#correct for multiple testing
sorted$adj <- p.adjust(as.numeric(sorted[,5]), "fdr")

fn_out="pval.count.table.txt_new"
write.table(sorted, file=paste(path,fn_out,sep=""),
            append=F, quote=F, 
            row.names=T, col.names=T, 
            dec=".", sep="\t")
