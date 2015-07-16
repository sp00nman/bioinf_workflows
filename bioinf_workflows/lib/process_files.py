import pandas as pd
import os
from utils import tools as ts


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
                     names=["chr", "start", "end",
                            "id", "mapq", "strand",
                            "cigar"])
    df['end'] = df['start']+1
    df.to_csv(outbedfile, sep="\t", index=0, header=0)

    return None


def count_insertions(exonfile,
                     intronfile,
                     outcountfile,
                     dn):
    """
    Count number of insertions within exons or introns
    :param exonfile: BED file with insertions within exons.
    :param intronfile: BED file with insertions within introns
    :param outcountfile: Plain text file with insertions counts for each gene.
    :param dn: path where Rscripts are located
    :return:Command to be executed; type str
    """

    cmd_count = "Rscript --vanilla " + dn \
                + "/" + "count_insertions.R %s %s %s" % (exonfile,
                                                             intronfile,
                                                             outcountfile)
    return cmd_count


def fisher_test(infile,
                control_file,
                outfile,
                dn):
    """
    Perform fisher-test for each gene between case and control screen.
    :param infile: input file
    :param control_file: control screen file
    :param outfile: output file
    :param dn: path where Rscripts are located
    :return: Command to be executed; type str
    """

    cmd_count = "Rscript --vanilla " + dn \
                + "/" + "fisher_test.R %s %s %s " % (infile,
                                                         control_file,
                                                         outfile)
    return cmd_count


def plot_results(infile,
                 refseq_file,
                 outfile,
                 dn):
    """
    Plot results in bubble plot with R.
    :param infile: input file
    :param refseq_file: refseq file
    :param outfile: output file
    :param dn: path where Rscripts are located
    :return: Command to be executed; type str
    """

    cmd_plot = "Rscript --vanilla " + dn \
               + "/" + "/plot_screen.R %s %s %s " % (refseq_file,
                                                        infile,
                                                        outfile)
    return cmd_plot


def browser_track(input_file,
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
                         names=["chr", "start", "end",
                                "id", "mapq", "strand",
                                "cigar", "ens_chr", "ens_start",
                                "ens_end", "ens_annotation",
                                "ens_strand", "ens_ensid",
                                "ens_gsymbol", "ens_transid",
                                "ens_num"])
        track_name = df['id'] + "~" + df['ens_gsymbol'] + "~" \
                     + annotation_name + df['strand'] + "/" + df['ens_strand']
        reshape = df.loc[:,['chr', 'start', 'end']]
        reshape['track_name'] = track_name
        reshape['score'] = 0
        reshape['strand'] = df['strand']

    else:
        df = pd.read_csv(input_file, sep="\t",
                         names=["chr", "start", "end",
                                "id", "mapq", "strand", "cigar",
                                "e_chr", "e_start", "e_end",
                                "e_peak_id", "e_score", "e_strand", "e_num"])
        reshape = df.loc[:, ['chr', 'start', 'end',
                             'e_peak_id', 'e_score', 'strand']]

    reshape.to_csv(output_file, mode='w', sep="\t", index=0, header=0)

    return None


def report_statistics():
    # number of mapped reads...duplicates,..insertion count..
    pass
