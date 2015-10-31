"""
Collection of functions to execute software for bamutils
"""


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