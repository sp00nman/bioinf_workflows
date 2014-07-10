Analysis workflow
=================


### Software requirements
 
+ [Picard] (http://picard.sourceforge.net/)
+ [Bowtie2] (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
+ [Bedtools] (https://github.com/arq5x/bedtools2)
+ [Samtools] (http://samtools.sourceforge.net/)


### Step-by-step workflow

Set environment variables. [Reference genome hg19 ](http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/)

```bash
export $REFERENCE="/path/to/reference"
export $EXONS="/path/to/exons"
export $INTRONS="/path/to/introns"
export $SCREEN="name_of_your_screen"

export $BOWTIE2="/path/to/bowtie"
export $NGS_PICARD="/path/to/picard"
export $INTERSECTBED="path/to/bedtools"
export $SAMTOOLS="path/to/samtools"

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

Filter for reads that

1. don't have a reported alignment (SAM flag = 4 )
2. have multiple alignments (XS:i)
3. have a mapping quality MAPQ higher than 20

SAM format and SAM flags explained.

+ [Description of SAM format](http://samtools.github.io/hts-specs/SAMv1.pdf)
+ [Samflags explained](http://picard.sourceforge.net/explain-flags.html)

```bash
samtools view -H -S ${SCREEN}.sam >${SCREEN}.filt.header.sam
awk -F"\t" '($0 ~ /^@/) {NR--; next}; ($2!=4 && ($13 !~/XS:i/) && $5>20)' ${SCREEN}.sam \
>> ${SCREEN}.filt.header.sam
```


Sort samfile by coordinate with picard.

```bash
java -jar $NGS_PICARD/SortSam.jar \
INPUT=${SCREEN}.filt.header.sam \
OUTPUT=${SCREEN}.filt.header.sorted.sam \
SORT_ORDER=coordinate
```

Mark & remove duplicates with picard.

```bash
java -jar $NGS_PICARD/MarkDuplicates.jar \
INPUT=${SCREEN}.filt.header.sorted.sam \
OUTPUT=${SCREEN}.filt.header.sorted.rem_dupl.sam \
METRICS_FILE=${SCREEN}.rem_dupl.metrics.txt \
REMOVE_DUPLICATES=true
```

Convert samfile to bedfile

```bash
awk 'BEGIN{OFS="\t"} \
!/^@/ \
{end=$4+1; print $3,$4,end,$1,$2,$5,$6,\
 $7,$8,$9,$10,$11,$12,$13,$14,$15,\
 $16,$17,$18,$19,$20 }' ${SCREEN}.filt.header.sorted.rem_dupl.sam \
 >${SCREEN}.filt.header.sorted.rem_dupl.bed
```

Remove insertions 1 or 2 base pairs away

```bash
sort -k1,1 -k2n ${SCREEN}.filt.header.sorted.rem_dupl.bed | 
awk '
BEGIN {OFS="\t"} 
{ROW[NR] = $2; sROW[NR] = $0}  \
END { for (i=1;i<=NR;i++) print sROW[i],ROW[i-1]} 
' |
awk '
BEGIN {OFS="\t"} 
{subtract = $2-$22; print $0, subtract}
' |
awk '!($23<=2 && $23>0)' |
cut -f1-21 >${SCREEN}.filt.header.sorted.rem_dupl.bp.bed
```

Annotate insertions. 

```bash
$INTERSECTBED \
-a ${SCREEN}.filt.header.sorted.rem_dupl.bp.bed \
-b ${EXONS} \
-wo >${SCREEN}.filt.header.sorted.rem_dupl.bp.exon.bed

$INTERSECTBED \
-a ${SCREEN}.filt.header.sorted.rem_dupl.bp.bed \
-b ${INTRONS} \
-wo >${SCREEN}.filt.header.sorted.rem_dupl.bp.intron.bed
```

Annotate insertions withing overlapping genes.

![overlapping](https://github.com/sp00nman/bionf_workflows/blob/master/img/overlapping.png?raw=true)

| Gene (strand)       | insertions (chr[n]:position |  strand  | intronic/exonic   | group       |
| :------------------ |:----------------------------|:---------|:------------------|:------------|
| GENE A (+)          | chr1:4,389,753              | +        | intronic          | mutagenic  |
| GENE B (-)          | chr1:4,389,753              | +        | intronic          | silent      |


```bash
awk '
(($5==16 && $26=="+") || \
($5==0 && $26=="-"))
' ${SCREEN}.filt.header.sorted.rem_dupl.bp.intron.bed \
>${SCREEN}.filt.header.sorted.rem_dupl.bp.intron.antisense.bed

awk '
(($5==16 && $26=="-") || \
($5==0 && $26=="+"))
' ${SCREEN}.filt.header.sorted.rem_dupl.bp.intron.bed \
>${SCREEN}.filt.header.sorted.rem_dupl.bp.intron.sense.bed

cat  ${SCREEN}.filt.header.sorted.rem_dupl.bp.exon.bed \
${SCREEN}.filt.header.sorted.rem_dupl.bp.intron.sense.bed \
>${SCREEN}.mutagenic.insertions.bed
```

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
| chr1:4,431,342             | -       | exonic            | mutagenci    |

(*) strand refers to how the insertions was mapped to the genome.

Count insertions.

```bash
cut -f28 ${SCREEN}.correct.insertions.bed | sort | uniq -c | \
sort -k1 -r -n | awk '{print $1"\t"$2}' >${SCREEN}.mutagenic.insertions.counts.bed

cut -f28 ${SCREEN}.filt.header.sorted.rem_dupl.bp.intron.antisense.bed | \
sort | uniq -c | sort -k1 -r -n | awk '{print $1"\t"$2}' \
>${SCREEN}.silent.insertions.counts.bed

awk -F"\t" '
NR==FNR {f1[$2]=$0; next}
($2 in f1) \
{print $0"\t"f1[$2]} \
' ${SCREEN}.mutagenic.insertions.counts.bed \
${SCREEN}.silent.insertions.counts.bed | \
awk -F"\t" '
BEGIN{OFS="\t"}
{ratio=$3/$1} {print $2,$1,$3,ratio}' >${SCREEN}.counts.table.txt

```

(Optional) Fisher test

| Gene A                 | SCREEN         | CONTROL |
| :--------------------- |:---------------|:--------|
| number of insertions   | 24             | 48      |
| no insertions          | 1303           | 1730    |


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


### Upload your data to ensemble 

![NF1 insertions](https://github.com/sp00nman/bionf_workflows/blob/master/img/NF1_2.png?raw=true)


### Summary table
[google spreadsheet with results of 41 screens](https://docs.google.com/spreadsheets/d/1XcimT1Aj45mjhUPsX4qHMlbFcj62h4InKtYPB6O1S4U/edit#gid=0)

### Plotting your data
[Plot your data using R](../doc/plot_data.md)


