"""
Collection of functions to execute software for picard
"""


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