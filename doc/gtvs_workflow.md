Analysis workflow
=================


### Software requirements
 
+ [Picard] (http://picard.sourceforge.net/)
+ [Bowtie2] (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
+ [Bedtools] (https://github.com/arq5x/bedtools2)


### Workflow

Set environment variables. [Reference genome hg19 ](http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/)

```bash
export $REFERENCE="/path/to/reference"
export $SCREEN="name_of_your_screen"
```

Align data to reference genome hg19. Parameter ```--sensitive``` equals to ```-D 15 -R 2 -L 22 -i S,1,1.15``` 

```bash
bowtie2 \
-p 2 \
--end-to-end \
--sensitive \
-x $REFERENCE \
-U ${SCREEN}.fastq \
-S ${SCREEN}.sam
```


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


### Summary table
[goole spreadsheet with results of 41 screens](https://docs.google.com/spreadsheets/d/1XcimT1Aj45mjhUPsX4qHMlbFcj62h4InKtYPB6O1S4U/edit#gid=0)

### Plotting your data


