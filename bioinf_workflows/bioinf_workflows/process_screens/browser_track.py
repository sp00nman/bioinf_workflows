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
    "e_ensids",
    "e_gsymbols",
    "e_annotation",
    "peak_length",
    "bp_overlap",
    "percent_overlap",
    "e_num"
]


def create_track(input_file,
                 output_file,
                 annotation_method):
    """
    Create browser tracks.
    :param input_file: input file
    :param output_file: output file
    :param annotation_name: name of annotation
    :return: None
    """

    if annotation_method == "exon-intron-bed":
        colnames = exon_intron_bed_BED_FILE
        track_genesymbol = 'ens_gsymbol'
        anno = 'ens_annotation'
        estrand = 'ens_strand'

    else:
        colnames = kbm7_bed_BED_FILE
        track_genesymbol = 'e_peak_id'
        anno = 'peak_length'
        estrand = 'e_num'

    df = pd.read_csv(input_file,
                     sep="\t",
                     names=colnames)

    track_name = df['id'] \
                 + "~" \
                 + df[track_genesymbol]\
                 + "~" \
                 + df[anno].apply(str) \
                 + "~" \
                 + df['strand'] \
                 + "/" + df[estrand].apply(str)

    reshape = df.loc[:, ['chr', 'start', 'end']]
    reshape['track_name'] = track_name
    reshape['score'] = 0
    reshape['strand'] = df['strand']

    out_filehandle = open(output_file, 'w')
    reshape.to_csv(out_filehandle,
                   sep="\t",
                   index=0,
                   header=0)
    out_filehandle.close()
    return None
