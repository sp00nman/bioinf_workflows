"""
Perform fisher-test
"""

import pandas as pd
import scipy.stats as spm
import numpy as np
import rpy2.robjects as R
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')

hg19_COUNT_COLUMN = [
    "gsymbol",
    "group_size",
    "total_size_minus_group_size",
    "ensids",
    "transids",
    "chr",
    "g_strand",
    "insertions_pos",
    "insertions_mapq",
    "insertions_strand",
    "insertions_annotation"
]

hg19_merge_COUNT_COLUMN = [
    'gsymbol_x',
    'group_size_x',
    'total_size_minus_group_size_x',
    'group_size_y',
    'total_size_minus_group_size_y',
    'pvalue',
    'fdr',
    'ensids',
    'transids_x',
    'chr_x',
    'g_strand_x',
    'insertions_pos_x',
    'insertions_mapq_x',
    'insertions_strand_x',
    'insertions_annotation_x'
]

kbm7_COUNT_COLUMN = [
    "gsymbol",
    "group_size",
    "total_size_minus_group_size",
    "peak_pos",
    "ensids",
    "gsymbols",
    "insertions_annotation",
    "peak_length",
    "bp_overlap",
    "percent_overlap",
    "insertions_pos"
]

kbm7_merge_COUNT_COLUMN = [
    'gsymbol',
    'group_size_x',
    'total_size_minus_group_size_x',
    'group_size_y',
    'total_size_minus_group_size_y',
    'pvalue',
    'fdr',
    "peak_pos_y",
    "ensids",
    "gsymbols_y",
    "insertions_annotation_y",
    "peak_length_y",
    "bp_overlap_y",
    "percent_overlap_y",
    "insertions_pos_y"
]


def fisher_test(group_size_data,
                total_size_minus_group_size_data,
                group_size_control,
                total_size_minus_group_size_control):

    oddsratio, pvalue = spm.fisher_exact(
        [[group_size_data, total_size_minus_group_size_data],
         [group_size_control, total_size_minus_group_size_control]])
    return pvalue


def read_files(filename, colnames):

    return pd.read_csv(filename,
                       sep="\t",
                       names=colnames,
                       dtype={'group_size': np.int32,
                              'total_size_minus_group_size': np.int32})


def analysis_workflow(infile,
                      control_file,
                      outfile,
                      annotation_method):

    if annotation_method == "exon-intron-bed":
        colnames = hg19_COUNT_COLUMN
        col_out = hg19_merge_COUNT_COLUMN
        merge_on = "ensids"
    else:
        colnames = kbm7_COUNT_COLUMN
        col_out = kbm7_merge_COUNT_COLUMN
        merge_on = "gsymbol"

    #read in data & control file
    data = read_files(infile, colnames)
    control = read_files(control_file, colnames)

    # merge table - inner == intersection
    merged_table = pd.merge(
        data, control,
        how="inner",
        on=merge_on)

    # drop duplicates
    merged_table = merged_table.drop_duplicates()

    # this is not necessary anymore, I guess
    # dropna_table = merged_table.dropna(how="any")

    # apply fisher-test for each row
    merged_table['pvalue'] = merged_table.apply(
        lambda row: fisher_test(
            row['group_size_x'],
            row['total_size_minus_group_size_x'],
            row['group_size_y'],
            row['total_size_minus_group_size_y']),
        axis=1
    )

    # fdr correction; for now only works with rpy2
    sorted_table = merged_table.sort('pvalue')
    sorted_table['fdr'] = stats.p_adjust(
        FloatVector(sorted_table['pvalue']),
        method='BH')

    out_handle = open(outfile, 'w')

    sorted_table.to_csv(
        out_handle,
        columns=col_out,
        sep="\t",
        index=0
    )

    out_handle.close()

    return None
