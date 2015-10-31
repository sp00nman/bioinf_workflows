"""
Perform fisher-test
"""

import pandas as pd
import scipy.stats as spm
import numpy as np
import rpy2.robjects as R

COUNT_COLUMN = [
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


def fisher_test(group_size_data,
                total_size_minus_group_size_data,
                group_size_control,
                total_size_minus_group_size_control):

    oddsratio, pvalue = spm.fisher_exact(
        [[group_size_data, total_size_minus_group_size_data],
         [group_size_control, total_size_minus_group_size_control]])
    return pvalue


def read_files(filename):

    return pd.read_csv(filename,
                       sep="\t",
                       names=COUNT_COLUMN,
                       dtype={'group_size': np.int32,
                              'total_size_minus_group_size': np.int32})


def analysis_workflow(infile,
                     control_file,
                     outfile):

    #read in data & control file
    data = read_files(infile)
    control = read_files(control_file)

    # merge table - inner == intersection
    merged_table = pd.merge(
        data, control,
        how="inner",
        on="ensids")

    # this is not necessary anymore, I guess
    # dropna_table = merged_table.dropna(how="any")

    # apply fisher test for each row
    merged_table['pvalue'] = merged_table.apply(
        lambda row: fisher_test(
            row['group_size_x'],
            row['total_size_minus_group_size_x'],
            row['group_size_y'],
            row['total_size_minus_group_size_y']),
        axis=1
    )

    # fdr correction; for now only works with rlibrary
    sorted_table = merged_table.sort('pvalue')
    # sorted_table['fdr'] = ...

    out_handle = open(outfile, 'w')
    sorted_table.to_csv(
        out_handle,
        sep="\t",
        index=0
    )
