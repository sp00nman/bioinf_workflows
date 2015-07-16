bionf_workflows
===============

Collection of bioinformatics workflows to analyze next-generation sequening data (Illumina Hiseq).

Installation procedure
======================

Clone the repository from github

```bash
git clone git@github.com:sp00nman/bioinf_workflows.git
cd bioinf_workflows
```

Install module:
```bash
python setup.py install 
```

or with `pip`:
```bash
pip install .
```

`--user` command can be added to install command if installation to home directory is necessary (eg. no root permission).

Python dependencies will be installed:
```bash
pandas
```


Analysis workflow
=================

## Software requirements
 
+ [Picard] (http://picard.sourceforge.net/)
+ [Bowtie2] (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
+ [Bedtools] (https://github.com/arq5x/bedtools2)
+ [Samtools] (http://samtools.sourceforge.net/)
+ [ngsutils - only for --stage statistics ](http://ngsutils.org/)

## Usage

```bash
usage: gt_screen_workflow.py [-h] [--debug DEBUG] [--stage STAGE]
                             [--project_name PROJECT_NAME]
                             [--output_dir OUTPUT_DIR]
                             [--sequences_dir SEQUENCES_DIR]
                             [--sample_file SAMPLE_FILE] [--genomes GENOMES]
                             [--genome_version GENOME_VERSION]
                             [--bowtie2 BOWTIE2] [--annotation ANNOTATION]
                             [--control_file CONTROL_FILE]
                             [--refseq_file REFSEQ_FILE] [--num_cpus NUM_CPUS]

Genetic screen workflow 0.0.1

optional arguments:
  -h, --help            show this help message and exit
  --debug DEBUG         Debug level
  --stage STAGE         Limit job submission to a particular analysis stage.[a
                        ll,(bam2fastq),alignment,filter,sort,duplicates,index,
                        insertions,annotate,count,fisher,plot,browser,statisti
                        cs]
  --project_name PROJECT_NAME
                        Name of project directory.
  --output_dir OUTPUT_DIR
                        Name of output directory.
  --sequences_dir SEQUENCES_DIR
                        Directory of sequence files.
  --sample_file SAMPLE_FILE
                        Input filename.
  --genomes GENOMES     Path to genome.
  --genome_version GENOME_VERSION
                        Genome version
  --bowtie2 BOWTIE2     Path to bowtie2 indices
  --annotation ANNOTATION
                        annotation file. Format (tab-separated: 
                        exon /path/to/annotation_1.txt 
                        intron /path/to/annotation_2.txt
  --control_file CONTROL_FILE
                        Control file with insertions for fisher-test.
  --refseq_file REFSEQ_FILE
                        Refseq file with start & end position of gene.
  --num_cpus NUM_CPUS   Number of cpus.

```

## --stage

### [all]
Process all analysis steps.

### [bam2fastq]
Convert unmapped bam files to fastq files with Picard (SamToFastq). Is not included if --stage == all.

### [alignment]

Align data to reference genome hg19 [Reference genome hg19 ](http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/). Parameter ```--sensitive``` equals to ```-D 15 -R 2 -L 22 -i S,1,1.15``` 

### [filter]

Filter for reads that

1. don't have a reported alignment (SAM flag = 4 )
2. have multiple alignments (XS:i)
3. have a mapping quality MAPQ higher than 20

SAM format and SAM flags explained.

+ [Description of SAM format](http://samtools.github.io/hts-specs/SAMv1.pdf)
+ [Samflags explained](http://picard.sourceforge.net/explain-flags.html)

### [sort]

Sort samfile by coordinate with [Picard] (http://picard.sourceforge.net/)

### [duplicates]

Mark & remove duplicates with [Picard] (http://picard.sourceforge.net/)

### [insertions]

Remove insertions 1 or 2 base pairs away. Insertions that are within close proximity could be a result of mapping errors due to sequencing errors and are therefore removed from the analysis.

### [annotate]

Annotate insertions. For annotation the genome assembly hg19 (ucsc) combined with the ensembl annotation (build 70) is used. The canonical transcripts for each gene is used as reference gene model. The bed file with coordinates for exons and introns is formated using the following [script] (https://gist.github.com/sp00nman/e9adb3c7e207c0de03d7)

Annotate insertions within overlapping genes.

![overlapping](https://github.com/sp00nman/bionf_workflows/blob/master/img/overlapping.png?raw=true 50x100)

| Gene (strand)       | insertions (chr[n]:position |  strand  | intronic/exonic   | group       |
| :------------------ |:----------------------------|:---------|:------------------|:------------|
| GENE A (+)          | chr1:4,389,753              | +        | intronic          | mutagenic  |
| GENE B (-)          | chr1:4,389,753              | +        | intronic          | silent      |

"Gene trap-based insertional mutagenesis operates by random insertion of a splice acceptor followed by a GFP
marker and termination sequence into the genome, thus disrupting gene expression. The insertion site must
be determined to identify the disrupted gene.....Of these, 321 clones contained gene-trap insertions
in a coding exon, directly disrupting the respective open reading frame. Insertions in introns are predictive
to be mutagenic if the gene-trap cassette is inserted in the sense orientation." T. Bürckstümmer et al., "a reversible gene trap collection empowers haploid genetics in human cells", Nature Methods, 2013.

![grouping](https://github.com/sp00nman/bionf_workflows/blob/master/img/grouping.png?raw=true 50x100)

| Gene A (chr[n]:position)   | strand(*)  | intronic/exonic   | group         |
| :------------------------- |:--------|:------------------|:--------------|
| chr1:4,389,753             | +       | intronic          | mutagenic    |
| chr1:4,399,100             | -       | intronic          | silent        |
| chr1:4,443,563             | +       | exonic            | mutagenic    |
| chr1:4,431,342             | -       | exonic            | mutagenic    |

(*) strand refers to how the insertions was mapped to the genome.

### [count]

Number of mutagenic insertions for each gene is counted. Equivalent to a UNIX [uniq -c] command.

### [fisher]

Each genetrap screen experiment represents an unreplicated data set and data analysis usually
proceeds on a gene-by-gene by organizing the data in a 2x2 table (table 1). Probably the most natural
test for significantely enriched inactivating insertions in a control vs experimental screen in
unreplicated data is Fisher's exact test.

Table 1

| Fisher-test                               | experiment   | control   | total    |
| :-----------------------------------------|:-------------|:----------|:---------|
| insertions (in gene A)                    | 11           | 76        | 87       |
| no insertions ([all other genes]-[gene A]) | 11,060       | 457,672   | 468,732  |
| total                                     | 11,071       | 457,748   | 468,819  |

Table 1 is an example of such a two-by-two or contingency table. For the given example the p-value is the probability of
observing a value 11 (11 insertions) or larger in gene A in the upper-left-hand cell in table 1,
assuming that the null hypothesis is true and that the four marginal totals 87, 468732, 11071,
457748 are given. Given a hypergeometric distribution the probability for gene A in table 1 is 6.75e-06. With a type I error of 5% the null hypothesis is rejected.
For each of the remaining 2528 genes a similar statistical test is performed. If for all of these 2528
genes the null hypothesis is actually true, 126 genes (2528*0.05) would turn up to be signicant
at a p-value below 0.05 just due to chance. These 126 genes are false positives and a procedure is
required to control the number of false-positives.

#### FDR-correction

Several procedure exists that attempt to assign an adjusted p-value to eacht test or similarly reduce
the p-value threshold. These adjustment methods include the Bonferroni correction ("Bonferroni")
in which the p-values are multiplied by the total number of comparisons (n). The procedure rejects
all hypotheses where the adjusted-pvalue is below a predened alpha. Less conservative corrections
were developed by Holm (1979) ("Holm"), Hochberg (1988) ("Hochberg"), Sidak (1967)("SidakSD")
and ("SidakSD"). These procedure are designed to give strong control of the family-wise error rate
(FWER). However, if n is very large, as is often the case in high-throughput screens including
genetrap screens, these procedures can lead to quite conservative results (insisting on a 5% error
rate).
For large n, procedures that control the false discovery rate (FDR) like Benjamini & Hochberg
(1995) ("BH"), and Benjamini & Yekutieli (2001) ("BY") are prefered. The FDR is a less stringent
condition than the FWER, consequently these methods are more powerful. Multiple testing
procedure were implemented in R using library multtest.
The function load.process.data calculates adjusted p-values for different multiple testing procedures
and assigns to each gene a chromosomal start and end position in the genome. Gene names
with chromosomal position of each gene is taken from a refseq annotation file.

It is pivotal to understand the concept of a false discovery
in contrast to a false positive. A false positive occurs when a null hypothesis (no significant
difference in number of insertions between control and case) is rejected when it is actually true.
Therefore the false positive rate is the proportion of false positives out of all true null hypotheses
that were incorrectly called significant.
A false discovery is also a false positive but as the name implies this time we are concerned with
the false positives among the set of comparisons that were called significant. When performing
many statistical tests the false discovery rate is more informative than the false positive rate. If
we are willing to accept a FDR of 5% this implies that among all test that were called significant,
about 5% will be false positives.

As an example, let's have a look at the previous output which shows the top 11 genes after multiple
testing correction odered by adjusted p-value ("BH"). There are 2528 genes in this experiment.
For gene X, the p-value ("rawp") is 0.000153 whereas the adjusted p-value (also
referred to as q-value) is 0.03856. The p-value of 0.000153 implies
a 0.0153% chance of false positives. Having 2528 genes we expect between 1 and 2 false positives
(2528*0.000153) on average. For this particular experiment, there are 11 genes with a p-value
smaller than 0.000153 or less, consequently 1 or 2 out of these will be false positives. Now let's
look at the q-value for BH which is a bit larger at 0.03856. This means, we would expect 3.856%
of all genes with a q-value less than 0.03856 to be false positives. This is very informative, now we
expect among all gene (11) with a q-value of less than 0.03856 (11*0.03856) 0,42416 false positives,
less than one false positive.
To put it in a nutshell, the false discovery rate gives a far better estimate of the level of false
positives for a threshold value.

### [plot]

Traditionally, genetrap screens were plotted as gene names (alphabetical order) on the x-axis vs
q-values on the y-axis. Using R library(ggplot2) genetrap screens can now be displayed with
their chromosomal position of each gene on the x-axis and the corresponding -log10 transformed
q-value on the log scale y-axis.

### [browser]

The output file with the following extension ${SCREEN}.mutagenic.genome_browser.bed can be visualized with a genome browser. 
+ [UCSC] genome browser (https://genome.ucsc.edu/cgi-bin/hgCustom?hgsid=382317369_ghDqdehOZl31AdcEk6KiErMEyhw8)
+ [Ensembl] genome browser (http://www.ensembl.org/)
![NF1 insertions](https://github.com/sp00nman/bionf_workflows/blob/master/img/NF1_2.png?raw=true 50x100)

### [statistics]

Not implemented yet.


### Output files explained

| file extension                        | explanation                                                           |
| :-------------------------------------| :---------------------------------------------------------------------|
| fastq                                 | FASTQ file                                                            |
| aligned.sam                           | [SAM file](http://samtools.github.io/hts-specs/SAMv1.pdf)             |
| aligned.bam                           | BAM file                                                              |
| filt_aligned.bam                      | BAM file, filtered according to [filter]                              |
| sorted_filt_aligned.bam               | BAM file, filtered and sorted                                         |
| reord_sorted_filt_aligned.bam         | BAM file, filtered, sorted and reordered                              |
| rmdupl_reord_sorted_filt_aligned.bam  | BAM file, after removal of duplicate reads                            |
| rmdupl_reord.sorted_filt_aligned.bai  | index of BAM file, after removal of duplicate reads                   |
| rmdupl_reord.sorted_filt_aligned.sam  | SAM file, after removal of duplicate reads                            |
| rm2bp_insertions.sam                  | SAM file, after removal of insertions 1 or 2 bp away (without header) |
| header.txt                            | header of BAM file                                                    |
| header_rm2bp_insertions.sam           | SAM file, after removal of insertions 1 or 2 bp away (with header)    |
| header_rm2bp_insertions.bam           | BAM file, after removal of insertions 1 or 2 bp away (with header)    |
| header_rm2bp_insertions.bed           | BED file                                                              |
| fix_header_rm2bp_insertions.bed       | BED file, added 1 bp to end position                                  |
| count_table.txt                       | Plain text file with counts of insertios for each gene                |
| fisher_test.txt                       | Plain text file with result table of fisher-test calculations         |
| bubble_plot.pdf                       | R generated bubble plot as pdf                                        |
| browser_track.txt                     | browser track files - can be uploaded to genome browser               |
| summary_statistics.txt                | not implemented yet                                                   |
| aligned.bai                           | index of aligned.bam file                                             |


###TODOs
+ gene lists of highly gene-trapped genes (?)
