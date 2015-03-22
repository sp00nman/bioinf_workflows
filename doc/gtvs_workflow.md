Analysis workflow
=================


## Software requirements
 
+ [Picard] (http://picard.sourceforge.net/)
+ [Bowtie2] (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
+ [Bedtools] (https://github.com/arq5x/bedtools2)
+ [Samtools] (http://samtools.sourceforge.net/)


## Usage

```bash
usage: gt_screen_workflow.py [-h] [--debug DEBUG] [--stage STAGE]
                             [--project_name PROJECT_NAME]
                             [--output_dir OUTPUT_DIR]
                             [--sequences_dir SEQUENCES_DIR]
                             [--sample_file SAMPLE_FILE] [--genomes GENOMES]
                             [--genome_version GENOME_VERSION]
                             [--bowtie2 BOWTIE2]
                             [--annotation_exon ANNOTATION_EXON]
                             [--annotation_intron ANNOTATION_INTRON]
                             [--control_file CONTROL_FILE]
                             [--refseq_file REFSEQ_FILE] [--num_cpus NUM_CPUS]

Genetic screen workflow 0.0.1

optional arguments:
  -h, --help            show this help message and exit
  --debug DEBUG         Debug level
  --stage STAGE         Limit job submission to a particular analysis stage.
                        [all,alignment,filter, sort, duplicates,
                        index, insertions, annotate, count, fisher, plot]
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
  --annotation_exon ANNOTATION_EXON
                        Exon annotation file.
  --annotation_intron ANNOTATION_INTRON
                        Intron annotation file
  --control_file CONTROL_FILE
                        Control file with insertions for fisher-test.
  --refseq_file REFSEQ_FILE
                        Refseq file with start & end position of gene.
  --num_cpus NUM_CPUS   Number of cpus.

```

## --stage

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

![overlapping](https://github.com/sp00nman/bionf_workflows/blob/master/img/overlapping.png?raw=true 100x200)

| Gene (strand)       | insertions (chr[n]:position |  strand  | intronic/exonic   | group       |
| :------------------ |:----------------------------|:---------|:------------------|:------------|
| GENE A (+)          | chr1:4,389,753              | +        | intronic          | mutagenic  |
| GENE B (-)          | chr1:4,389,753              | +        | intronic          | silent      |

"Gene trap-based insertional mutagenesis operates by random insertion of a splice acceptor followed by a GFP
marker and termination sequence into the genome, thus disrupting gene expression. The insertion site must
be determined to identify the disrupted gene.....Of these, 321 clones contained gene-trap insertions
in a coding exon, directly disrupting the respective open reading frame. Insertions in introns are predictive
to be mutagenic if the gene-trap cassette is inserted in the sense orientation." T. Bürckstümmer et al., "a reversible gene trap collection empowers haploid genetics in human cells", Nature Methods, 2013.

![grouping](https://github.com/sp00nman/bionf_workflows/blob/master/img/grouping.png?raw=true )

| Gene A (chr[n]:position)   | strand(*)  | intronic/exonic   | group         |
| :------------------------- |:--------|:------------------|:--------------|
| chr1:4,389,753             | +       | intronic          | mutagenic    |
| chr1:4,399,100             | -       | intronic          | silent        |
| chr1:4,443,563             | +       | exonic            | mutagenic    |
| chr1:4,431,342             | -       | exonic            | mutagenic    |

(*) strand refers to how the insertions was mapped to the genome.

### [count]

Count insertions.

### [fisher]

### [plot]

### Output files explained

#### ${SCREEN}.filt.header.sorted.sam
+ [Description of SAM format](http://samtools.github.io/hts-specs/SAMv1.pdf)

##### ${SCREEN}.filt.header.sorted.rem_dupl.bp.exon.bed
|desription|example|
|:----------|:-----|
|chromosome (where insertions mapped to)|chr1|    
|start position (of insertion)|713901|  
|end position (of insertions)|713902 | 
|[QNAME](http://samtools.github.io/hts-specs/SAMv1.pdf)|HWI-ST815:70:D1F2MACXX:2:2303:6014:198907|       
|[FLAG](http://samtools.github.io/hts-specs/SAMv1.pdf)|16 |      
|[MAPQ](http://samtools.github.io/hts-specs/SAMv1.pdf)|40 |     
|[CIGAR](http://samtools.github.io/hts-specs/SAMv1.pdf)|51M |    
|[RNEXT](http://samtools.github.io/hts-specs/SAMv1.pdf)|*  |     
|[PNEXT](http://samtools.github.io/hts-specs/SAMv1.pdf)|0  |     
|[TLEN](http://samtools.github.io/hts-specs/SAMv1.pdf)|0  |     
|[SEQ](http://samtools.github.io/hts-specs/SAMv1.pdf)|CGGCAACCCACAGGTCCTGGCGGGGACGTCACTCTTACCAGTCCCCACTCT | 
|[QUAL](http://samtools.github.io/hts-specs/SAMv1.pdf)|###################A:CF?1:)AFA<<@BF@FA=HFCDDBB=?@?? |    
|[OPTIONAL FIELDS](http://samtools.github.io/hts-specs/SAMv1.pdf)|MD:Z:7A5A0A11G24 PG:Z:MarkDuplicates XG:i:0  NM:i:4  XM:i:4  XN:i:0  XO:i:0  AS:i:-8 YT:Z:UU | 
|chromosome of annotated gene |chr1  |  
|start position of feature|713664 | 
|end position of feature |714006 | 
|description of feature|exon  |  
|strand direction of feature|-  |     
|ensembl id associated to that feature|ENSG00000228327 |
|gene symbol associated to that feature|RP11-206L10.2  | 
|transcript id associated to that feature|ENST00000428504 |
|depricated, will be removed in later versions |1 |

### Upload your data to ensemble 
The output file with the following extension ${SCREEN}.mutagenic.genome_browser.bed can be visualized with a genome browser. 
+ [UCSC] genome browser (https://genome.ucsc.edu/cgi-bin/hgCustom?hgsid=382317369_ghDqdehOZl31AdcEk6KiErMEyhw8)
+ [Ensembl] genome browser (http://www.ensembl.org/)
![NF1 insertions](https://github.com/sp00nman/bionf_workflows/blob/master/img/NF1_2.png?raw=true)

### Summary table
[google spreadsheet with results of 41 screens](https://docs.google.com/spreadsheets/d/1XcimT1Aj45mjhUPsX4qHMlbFcj62h4InKtYPB6O1S4U/edit#gid=0)

### Parameter description
| metrics                 | description          |
| :----------------------- |:----------------------|
| TOTAL_NUMBER_READS      | total amount of sequencing reads |
| UNMAPPED_READS          | number of reads that map 0 times     |
| UNMAPPED_READS_PCT      | percentage of unmapped reads    |
| UNIQUELY_MAPPED_READS   | number of reads that map exactely 1 time |
| FILTERED_READS          | number of reads after filtering | 
| DUPLICATION_RATE_PCT    | percentage of duplicated reads |
| UNIQUE_READS            | number of reads after removal of duplicates |
| TOTAL_NUMBER_INSERTIONS | number of reads after removal of reads 1 or 2 base pairs away | 
| NUMBER_TARGETED_GENES   | number of gene that are targeted |
| NUMBER_INSERTIONS_INTERGENIC | number of insertions that are in between genes | 
| NUMBER_INSERTIONS_INTRAGENIC | number of insertions within introns and exons | 
| NUMBER_INSERTIONS_INTRONIC_SENSE | number of insertions that are within introns and sense | 
| NUMBER_INSERTIONS_INTRONIC_ANTISENSE | number of insertions that are within introns and antisense |
| NUMBER_INSERTIONS_EXONIC | number of insertions that are within exons |
| TOP10_GENES | number of top 10 genes |

### Some example results

![TOTAL_NUMBER_READS vs UNIQUE_READS](https://github.com/sp00nman/bionf_workflows/blob/master/img/img2.png?raw=true "TOTAL_NUMBER_READS vs UNIQUE_READS")

![DUPLICATION_RATE vs UNIQUE_READS](https://github.com/sp00nman/bionf_workflows/blob/master/img/img1.png?raw=true "DUPLICATION_RATE vs UNIQUE_READS")


### Visualize your data
+ [Plot bubble plot using R](../doc/plot_data.md)
+ [Create publication ready circos plots](../doc/circos_plots.md)


###TODOs
+ automate bubble plots
+ implement bash script to a python workflow
+ script for submitting jobs to SLURM
+ gene lists of highly gene-trapped genes (?)
