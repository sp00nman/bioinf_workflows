"""
Collection of functions to execute software for samtools
"""


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