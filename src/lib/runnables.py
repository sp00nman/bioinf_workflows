"""
Collection of functions that execute external software or UNIX commands.
"""


def bam2fastq(bamfile,
              fastqfile):
    """
    Converts BAM formatted files to FASTQ formatted files with PICARD.
    :param bamfile: Sample filename
    :param fastqfile:project sub-directory name
    :return: Command to be executed; type str
    """
    cmd_bam2fastq = "java -Xmx6g -jar $NGS_PICARD/SamToFastq.jar " \
                    "INPUT=%s " \
                    "FASTQ=%s" % (bamfile, fastqfile)
    return cmd_bam2fastq


def alignment(genome_path,
              fastqfile,
              samfile,
              num_cpus,
              metrics):
    """
    Alignment of reads to the reference genome with bowtie2.
    :param genome_path: version of the genome, eg [hg19, mm10]
    :param fastqfile: name of FASTQ formatted file
    :param samfile: name of output file
    :param num_cpus: for multi-threading purposes (-p option; bowtie2)
    :return: Command to be executed; type str
    """
    cmd_align = "bowtie2 " \
                "-p %s " \
                "--end-to-end " \
                "--sensitive " \
                "-x %s " \
                "-U %s " \
                "-S %s " \
                "--met-file %s" \
                % (num_cpus,
                   genome_path,
                   fastqfile,
                   samfile,
                   metrics)
    return cmd_align


def sam2bam(insamfile,
            outbamfile):
    """
    Convert SAM formatted file to BAM formatted file.
    :param insamfile: Input samfile
    :param outbamfile: Output bamfile
    :return: Command to be executed; type str
    """
    cmd_sam2bam = "samtools view " \
                  "-S " \
                  "-b %s " \
                  ">%s" % (insamfile,
                           outbamfile)
    return cmd_sam2bam


def filter_reads(inbamfile,
                 outbamfile,
                 mapq,
                 flag):
    """
    Filter read based on
    -q mapping quality of read
    -F only uniquely mapped reads
    :param inbamfile: name of BAM formatted file
    :param outbamfile: name of output file (BAM formatted and filtered)
    :param mapq: mapping quality
    :return: Command to be executed; type str
    """
    cmd_filter = "samtools view " \
                 "-b " \
                 "-q %s " \
                 "-F %s " \
                 "%s >%s" % (mapq,
                             flag,
                             inbamfile,
                             outbamfile)
    return cmd_filter


def sort_bam(inbamfile,
             outbamfile,
             sort_order="coordinate"):
    """
    Sort BAM file by variable
    :param inbamfile: name of BAM formatted file
    :param outbamfile: name of output file (BAM formatted and sorted)
    :param sort_order: sort by (default: coordinate)
    :return: Command to be executed; type str
    """

    cmd_sort = "java -Xmx6g -jar $NGS_PICARD/SortSam.jar " \
               "INPUT=%s " \
               "OUTPUT=%s " \
               "SORT_ORDER=%s" % (inbamfile,
                                  outbamfile,
                                  sort_order)
    return cmd_sort


def reorder_sam(inbamfile,
                outbamfile,
                genome_path):
    """
    Reorder BAM file.
    :param inbamfile: name of BAM formatted file
    :param outbamfile: name of output file
    :param genomes: path to genome
    :return: Command to be executed; type str
    """

    cmd_reorder = "java -Xmx6g -jar $NGS_PICARD/ReorderSam.jar " \
                  "INPUT=%s " \
                  "OUTPUT=%s " \
                  "REFERENCE=%s.fa" % (inbamfile,
                                    outbamfile,
                                    genome_path)
    return cmd_reorder


def count_duplicates(inbamfile,
                     out_statistics):
    """
    Output of the number of reads at each position. This is actually
    the number of duplicate reads at each position. If a position has
    multiple reads mapped to it, but they are not pcr duplicates, then
    there each will be reported separately.
    :param inbamfile: name of BAM formatted file
    :param out_statistics: name of output file
    :return: Command to be executed; type str
    """

    cmd_countdup = "bamutils pcrdup " \
                   + "-frag " \
                   + "-counts %s %s " % (out_statistics,
                                         inbamfile)
    return cmd_countdup


def remove_duplicates(inbamfile,
                      outbamfile,
                      metrics_file):
    """
    Remove duplicate reads.
    :param inbamfile: name of BAM formatted file
    :param outbamfile: name of output file
    :param metrics_file: file with summary statistics about
    :return: Command to be executed; type str
    """

    cmd_rmdup = "java -Xmx6g -jar $NGS_PICARD/MarkDuplicates.jar " \
                "INPUT=%s " \
                "OUTPUT=%s " \
                "METRICS_FILE=%s " \
                "REMOVE_DUPLICATES=true" % (inbamfile,
                                            outbamfile,
                                            metrics_file)
    return cmd_rmdup


def index_bam(inbamfile):
    """
    Index BAM file.
    :param inbamfile: name of BAM formatted file
    :return: Command to be executed; type str
    """

    cmd_index_bam = "samtools index %s" % (inbamfile)
    return cmd_index_bam


def bam2sam(inbamfile,
            outsamfile):
    """
    Convert BAM file to SAM file.
    :param inbamfile: name of BAM formatted file
    :param outsamfile: name of SAM file
    :return: Command to be executed; type str
    """

    cmd_bam2sam = "samtools view %s > %s" % (inbamfile,
                                             outsamfile)
    return cmd_bam2sam


def get_header(inbamfile,
               header):
    """
    Extract header from BAM file.
    :param inbamfile: name of BAM formatted file
    :param header: header file
    :return: Command to be executed; type str
    """

    cmd_get_header = "samtools view -H %s > %s" % (inbamfile,
                                                   header)
    return cmd_get_header


def concatenate_files(file_1,
                      file_2,
                      output_file):
    """
    Concatenate two files write output to new file. Uses standard UNIX cat
    :param file_1: 1st parameter
    :param file_2: 2nd parameter
    :param output_file: name of output file
    :return: Command to be executed; type str
    """

    cmd_concatenate_files = "cat %s %s >%s" % (file_1,
                                               file_2,
                                               output_file)
    return cmd_concatenate_files


def sam2bed(insamfile,
            outbedfile):
    """
    Convert SAM file to BED file
    :param insamfile: input SAM file
    :param outbedfile: input BED file
    :return: Command to be executed; type str
    """

    cmd_sam2bed = "bamToBed -cigar -i %s >%s " % (insamfile,
                                                  outbedfile)
    return cmd_sam2bed


def intersectbed(inbedfile,
                 annotation_file,
                 outbedfile):
    """
    Intersect two BED files by coordinate.
    :param inbedfile: input BED file
    :param annotation_file: annotation file
    :param outbedfile: output BED file
    :return: Command to be executed; type str
    """

    cmd_intersect = "intersectBed " \
                    "-a %s " \
                    "-b %s " \
                    "-wo >%s.bed" % (inbedfile,
                                     annotation_file,
                                     outbedfile)
    return cmd_intersect
