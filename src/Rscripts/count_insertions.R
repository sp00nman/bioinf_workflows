## count number of insertions within exons or introns

read.count.tables <- function(file){
  # read in count table
  header = c("chr", "start", "end",
             "id", "mapq", "strand", 
             "cigar", "chr_feature", "start_feature",
             "end_feature", "feature", "strand_feature", 
             "ensembl_gene", "genesymbol", "ensembl_transcript", 
             "num")
  df <- read.table(file, sep="\t")
  colnames(df) <- header
  return(df)  
}

separate.silent.mutagenic.insertions <- function(df_insertions){
  # diff orientation --> muatgenic
  # same orientation --> silent
  mutagenic <- df_insertions[df_insertions$strand == df_insertions$strand_feature,]
  silent <- df_insertions[df_insertions$strand != df_insertions$strand_feature,]
  return(list(mutagenic, silent))
}

count.insertions <- function(count_exon, count_intron){
  # merge
  # count only unique insertions
  # 
}


options <- commandArgs(trailingOnly=TRUE)
print options
#read.count.tables(options[1], options[2])
exon_f = "../../test/TEST.insertions.exon.bed"
intron_f = "../../test/TEST.insertions.intron.bed"

intron_insertions <- read.count.tables(intron_f)
exon_insertions <- read.count.tables(exon_f)

orientation <- separate.silent.mutagenic.insertions(intron_insertions)

