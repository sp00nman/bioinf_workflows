"""
Count number of insertions for each ensid (!!!)
Careful, I found cases with same genesymbol but different ensids....
"""

import pandas as pd

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


def count_total_uniq_insertions(data, columnid):
    """
    Count unique number of insertions based on illumina sequence identifier.
    :param data: pandas dataframe
    :param columnid: column id to use
    :return: number with total counts
    """
    return len(data[columnid].unique())


def count_insertions_pythonic(inbedfile,
                              outfile,
                              annotation_method):

    out_handle = open(outfile, 'a')

    if annotation_method == 'exon-intron-bed':
        dta = pd.read_csv(
            inbedfile,
            sep="\t",
            header=False,
            names=BED_FILE_ANNO)

        # calculate total unique insertions
        sum_uniq_ins = count_total_uniq_insertions(dta, columnid='id')

        # only select mutagenic insertions
        # all exons, and introns only if in sense orientation
        gb_mutagenic = dta[(dta.ens_annotation == "exon") |
                           ((dta.ens_annotation == "intron") &
                            (dta.strand == dta.ens_strand))]

        # groups all genes by
        gb = gb_mutagenic.groupby(['ens_ensid',
                                   'ens_gsymbol',
                                   'ens_transid'],
                                  sort=False,
                                  as_index=False)

        for (t1, t2, t3), group in gb:
            # colnames: ens_gsymbol,group_size,sum_uniq:ins-group_size,
            # ensid,transid,chr,ins_start,ins_mapq,ins_strand,annotation
            # TODO: implement this
            # this might make it easier
            # handle = open("asd.csv", "w")
            # wr = csv.writer(handle, delimiter="\t")
            # l is a list of lists
            # wr.writerows(l)
            # handle.close
            # #','.join([str(x) for x in list_num])

            line = t2 + "\t" \
                   + str(len(group)) + "\t" \
                   + str(sum_uniq_ins - len(group)) + "\t" \
                   + t1 + "\t" \
                   + t3 + "\t" \
                   + str(group['chr'].tolist()[0]) + "\t" \
                   + str(group['ens_strand'].tolist()[0]) + "\t" \
                   + str(group['start'].tolist()) + "\t" \
                   + str(group['mapq'].tolist()) + "\t" \
                   + str(group['strand'].tolist()) + "\t" \
                   + str(group['ens_annotation'].tolist()) + "\n"

            # print line to file and perform fisher test
            out_handle.write(line)

    # else annotation_method == 'locus-bed'
    else:
        dta = pd.read_csv(
            inbedfile,
            sep="\t",
            header=False,
            names=kbm7_bed_BED_FILE)

        # calculate total unique insertions
        sum_uniq_ins = count_total_uniq_insertions(dta, columnid='id')

        # groups all genes by peak-id
        gb = dta.groupby(['e_peak_id'],
                         sort=False,
                         as_index=False)

        for name, group in gb:
            line = name + "\t" \
                  + str(len(group)) + "\t" \
                  + str(sum_uniq_ins - len(group)) + "\t" \
                  + str(group['chr'].tolist()[0]) + ":" \
                  + str(group['e_start'].tolist()[0]) + "-" \
                  + str(group['e_end'].tolist()[0]) + "\t" \
                  + str(group['e_ensids'].tolist()[0]) + "\t" \
                  + str(group['e_gsymbols'].tolist()[0]) + "\t" \
                  + str(group['e_annotation'].tolist()[0]) + "\t"\
                  + str(group['peak_length'].tolist()[0]) + "\t" \
                  + str(group['bp_overlap'].tolist()[0]) + "\t" \
                  + str(group['percent_overlap'].tolist()[0]) + "\t" \
                  + ','.join([str(x) for x in group['start'].tolist()]) + "\n"

            out_handle.write(line)

    out_handle.close()

    return None

# some useful pieces of code
#for name, group in gb:
#    print str(group['chr'].tolist()[0]) + "\t" \
#          + str(group['start'].tolist()[0]) + "\t" \
#          + str(group['end'].tolist()[0]) + "\t" \
#          + name + "\t" + ','.join(group['e_ensgene'].tolist()) + "\t" \
#          + ','.join(group['e_gsymbol'].tolist()) + "\t" \
#          + ','.join(group['e_annotation'].tolist()) + "\t" \
#          + str(group['peak_length'].tolist()[0]) + "\t" \
#          + ','.join([str(x) for x in group['bp_overlap'].tolist()]) + "\t" \
#          + ','.join([str(x) for x in group['percent_overlap'].tolist()]) + "\n"
