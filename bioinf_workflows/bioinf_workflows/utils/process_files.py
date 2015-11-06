import pandas as pd

from bioinf_workflows.utils import tools as ts

BED_FILE = [
    "chr",
    "start",
    "end",
    "id",
    "mapq",
    "strand",
    "cigar"
]

BED_FILE_ANNO = [
    "chr",
    "start",
    "end",
    "id",
    "mapq",
    "strand",
    "cigar",
    "ens_chr",
    "ens_start",
    "ens_end",
    "ens_annotation",
    "ens_strand",
    "ens_ensid",
    "ens_gsymbol",
    "ens_transid",
    "ens_num"
]


def remove2bpinsertions(insamfile,
                        sam_out):
    """
    Remove insertions that are one or two base pairs away from each other.
    :param insamfile: input SAM file
    :param sam_out: file object to write content
    :return: file object to write content (???)
    """

    sam_file = ts.load_file(insamfile)

    for index in range(len(sam_file)):
        if index == 0:
            sam_out.writelines('\t'.join(sam_file[index-1]))
            sam_out.writelines("\n")
            continue

        # samfile chromosome position is column 4
        position = int(sam_file[index-1][3])
        next_read = sam_file[index]
        next_position = int(next_read[3])
        bp_distance = abs(next_position-position)

        if not (bp_distance <= 2):
            sam_out.writelines('\t'.join(sam_file[index]))
            sam_out.writelines("\n")

    return sam_out


def fix_end_position(inbedfile,
                     outbedfile):
    """
    Add 1 base pair to end position in BED formatted file
    :param inbedfile: input BED file
    :param outbedfile: output BED file
    :return: None
    """

    df = pd.read_csv(inbedfile, sep="\t",
                     names=BED_FILE)
    df['end'] = df['start']+1
    df.to_csv(outbedfile, sep="\t", index=0, header=0)

    return None




