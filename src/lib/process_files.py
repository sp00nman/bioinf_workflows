import pandas as pd
import os
from utils.tools import load_files


def remove2bpinsertions(project_name,
                        project_dir,
                        sample_file,
                        file_ext):

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "." + file_ext
    msg_rm2bpins = "Remove insertions 1 or 2 bp away."

    sam_file = load_files(input_file)
    sam_out = open(output_file, 'w')

    for index in range(len(sam_file)):
        #print index
        if index == 0:
            #print '\t'.join(sam_file[index])
            sam_out.writelines('\t'.join(sam_file[index-1]))
            sam_out.writelines("\n")
            continue

        # samfile chromosome position is column 4
        position = int(sam_file[index-1][3])

        next_read = sam_file[index]
        next_position = int(next_read[3])

        # bpdis --> base pair distance
        bp_dist = abs(next_position-position)

        if not (bp_dist <= 2):
            #print '\t'.join(sam_file[index])
            sam_out.writelines('\t'.join(sam_file[index]))
            sam_out.writelines("\n")

    sam_out.close()
    return msg_rm2bpins


def fix_end_position(project_name,
                     project_dir,
                     sample_file,
                     file_ext):

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "." + file_ext

    msg_fix = "Fix end position to be start+1."
    df = pd.read_csv(input_file, sep="\t",
                     names=["chr", "start", "end",
                            "id", "mapq", "strand", "cigar"])
    df['end'] = df['start']+1
    df.to_csv(output_file, sep="\t", index=0, header=0)

    return msg_fix


def count_insertions(project_name,
                     project_dir,
                     file_ext):

    exon = project_dir + "/" + project_name + "." + "insertions" + "." \
           + "exon" + "." + "bed"
    intron = project_dir + "/" + project_name + "." + "insertions" + "." \
             + "intron" + "." + "bed"
    output_file = project_dir + "/" + project_name + "." + file_ext

    dn = os.path.dirname(os.path.realpath(__file__))
    msg_count = "Count number of insertions."
    cmd_count = "Rscript --vanilla " + dn \
                + "/Rscripts/count_insertions.R %s %s %s" % (exon, intron,
                                                             output_file)
    return msg_count, cmd_count


def fisher_test(project_name,
                project_dir,
                sample_file,
                control_file,
                file_ext):

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "." + file_ext
    dn = os.path.dirname(os.path.realpath(__file__))
    msg_count = "Fisher-test for differential number of insertions."
    cmd_count = "Rscript --vanilla " + dn \
                + "/Rscripts/fisher_test.R %s %s %s " % (input_file,
                                                         control_file,
                                                         output_file)
    return msg_count, cmd_count


def plot_results(project_name,
                 project_dir,
                 refseq_file,
                 sample_file,
                 file_ext):

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "." + file_ext
    dn = os.path.dirname(os.path.realpath(__file__))
    msg_plot = "Create bubble plot for genetic screen results. "
    cmd_plot = "Rscript --vanilla " + dn \
               + "/Rscripts/plot_screen.R %s %s %s " % (refseq_file,
                                                        input_file,
                                                        output_file)
    return msg_plot, cmd_plot


def browser_track(project_name,
                  project_dir,
                  annotation_name,
                  file_ext):

    input_file = project_dir + "/" + project_name + "." + "insertions" + "." \
        + annotation_name + ".bed"
    output_file = project_dir + "/" + project_name + "." \
        + annotation_name + "." + file_ext

    msg_track = "Create browser tracks. "

    if annotation_name=="exon" or annotation_name=="intron":
        df = pd.read_csv(input_file, sep="\t",
                         names=["chr", "start", "end",
                                "id", "mapq", "strand",
                                "cigar", "ens_chr", "ens_start",
                                "ens_end", "ens_annotation",
                                "ens_strand", "ens_ensid",
                                "ens_gsymbol", "ens_transid",
                                "ens_num"])
        track_name = df['id'] + "~" + df['ens_gsymbol'] + "~" + annotation_name \
                     + df['strand'] + "/" + df['ens_strand']
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
        reshape = df.loc[:,['chr', 'start', 'end','e_peak_id','e_score','strand']]

    reshape.to_csv(output_file, mode='w', sep="\t", index=0, header=0)

    return msg_track


def report_statistics():
    # number of mapped reads...duplicates,..insertion count..
    pass