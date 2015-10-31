"""
Write browser track to upload to ucsc or ensembl.
"""

import pandas as pd

exon_intron_bed_BED_FILE = [
    "chr",
    "start",
    "end",
    "id",
    "mapq",
    "strand",
    "cigar"
]

kbm7_bed_BED_FILE = [
    "chr",
    "start",
    "end",
    "id",
    "mapq",
    "strand",
    "cigar",
    "e_chr",
    "e_start",
    "e_end",
    "e_peak_id",
    "e_score",
    "e_strand",
    "e_num"
]


def create_track(input_file,
                 output_file,
                 annotation_name):
    """
    Create browser tracks.
    :param input_file: input file
    :param output_file: output file
    :param annotation_name: name of annotation
    :return: None
    """

    if annotation_name == "exon" or annotation_name == "intron":
        df = pd.read_csv(input_file, sep="\t",
                         names=exon_intron_bed_BED_FILE)
        track_name = df['id'] + "~" + df['ens_gsymbol'] + "~" \
                     + annotation_name + df['strand'] + "/" + df['ens_strand']
        reshape = df.loc[:,['chr', 'start', 'end']]
        reshape['track_name'] = track_name
        reshape['score'] = 0
        reshape['strand'] = df['strand']

    else:
        df = pd.read_csv(input_file, sep="\t",
                         names=kbm7_bed_BED_FILE)
        reshape = df.loc[:, ['chr', 'start', 'end',
                             'e_peak_id', 'e_score',
                             'strand']]

    reshape.to_csv(output_file, mode='w', sep="\t", index=0, header=0)

    return None
