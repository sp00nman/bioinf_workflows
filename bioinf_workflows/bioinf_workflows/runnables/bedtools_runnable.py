"""
Collection of functions to execute software for bedtools
"""


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
                    "-wo >%s" % (inbedfile,
                                     annotation_file,
                                     outbedfile)
    return cmd_intersect
