Analysis workflow
=================


### Software requirements
 
+ [Picard] (http://picard.sourceforge.net/)
+ [Bowtie2] (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
+ [Bedtools] (https://github.com/arq5x/bedtools2)
+ [Samtools] (http://samtools.sourceforge.net/)


### Workflow

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

1. don't have a reported alignment (column 2 equals for 4 )
2. have multiple alignments (XS:i)
3. have a mapping quality MAPQ higher than 20

What are FLAGS? 
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
intersectBed \
-a ${SCREEN}.filt.header.sorted.rem_dupl.bp.bed \
-b ${EXONS} \
-wo >${SCREEN}.filt.header.sorted.rem_dupl.bp.exon.bed

 intersectBed \
-a ${SCREEN}.filt.header.sorted.rem_dupl.bp.bed \
-b ${INTRONS} \
-wo >${SCREEN}.filt.header.sorted.rem_dupl.bp.intron.bed

```

Assign insetions to groups

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
>${SCREEN}.correct.insertions.bed
```

Count insertions.

```bash
cut -f28 ${SCREEN}.correct.insertions.bed | sort | uniq -c | \
sort -k1 -r -n | awk '{print $1"\t"$2}' >${SCREEN}.correct.insertions.counts.bed

cut -f28 ${SCREEN}.filt.header.sorted.rem_dupl.bp.intron.antisense.bed | \
sort | uniq -c | sort -k1 -r -n | awk '{print $1"\t"$2}' \
>${SCREEN}.incorrect.insertions.counts.bed

awk -F"\t" '
NR==FNR {f1[$2]=$0; next}
($2 in f1) \
{print $0"\t"f1[$2]} \
' ${SCREEN}.correct.insertions.counts.bed \
${SCREEN}.incorrect.insertions.counts.bed | \
awk -F"\t" '
BEGIN{OFS="\t"}
{ratio=$3/$1} {print $2,$1,$3,ratio}' >${SCREEN}.counts.table.txt

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
[google spreadsheet with results of 41 screens](https://docs.google.com/spreadsheets/d/1XcimT1Aj45mjhUPsX4qHMlbFcj62h4InKtYPB6O1S4U/edit#gid=0)

### Plotting your data


