#!/usr/bin/env python

import argparse
import re
import logging
import os

from bioinf_workflows.runnables import bamutils_runnables as bamutils
from bioinf_workflows.runnables import bedtools_runnable as bedtools
from bioinf_workflows.runnables import samtools_runnables as samtools
from bioinf_workflows.runnables import picard_runnables as picard
from bioinf_workflows.runnables import bowtie_runnables as bowtie
from bioinf_workflows.runnables import R_runnables as rexe
from bioinf_workflows.process_screens import count_insertions
from bioinf_workflows.process_screens import fisher_test
from bioinf_workflows.process_screens import browser_track
from bioinf_workflows.utils import process_files as pf
from bioinf_workflows.utils import tools as ts


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Genetic screen workflow 0.0.1')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False,
                        help='Limit job submission to a particular '
                             'analysis stage.'
                             '[all,(bam2fastq),alignment,filter,sort,duplicates,'
                             'index,insertions,annotate,count,fisher,plot,'
                             'browser,statistics]')

    # directory and sample input
    parser.add_argument('--project_name', dest='project_name', required=False,
                        help='Name of project directory.')
    parser.add_argument('--output_dir', dest='output_dir', required=False,
                        help='Name of output directory.')
    parser.add_argument('--sequences_dir', dest='sequences_dir', required=False,
                        help='Directory of sequence files.')
    parser.add_argument('--sample_file', dest='sample_file', required=False,
                        help='Input filename.')

    # genome path & annotation
    parser.add_argument('--genomes', dest='genomes', required=False,
                        help='Path to genome.')
    parser.add_argument('--genome_version', dest='genome_version', required=False,
                        help='Genome version')
    parser.add_argument('--bowtie2', dest='bowtie2', required=False,
                        help='Path to bowtie2 indices')
    parser.add_argument('--annotation', dest='annotation', required=False,
                        help='annotation file. Format (tab-separated: '
                             'exon /path/to/annotation_1.txt '
                             'intron /path/to/annotation_2.txt)')
    parser.add_argument('--control_file', dest='control_file', required=False,
                        help='Control file with insertions for fisher-test.')

    # resources & parameters for plotting
    parser.add_argument('--refseq_file', dest='refseq_file', required=False,
                        help='Refseq file with start & end position of gene.')
    parser.add_argument('--ins_annotation', dest='ins_annotation',
                        required=False, help='gtf file for insertion annotation')
    parser.add_argument('--plot_option', dest='plot_option', required=False,
                        help='[png|pdf]')
    parser.add_argument('--fdr_cutoff', dest='fdr_cutoff', required=False,
                        help='FDR cutoff for plotting')

    # hardware related requirements
    parser.add_argument('--num_cpus', dest='num_cpus', required=False,
                        help='Number of cpus.')

    args = parser.parse_args()

    # set defaults
    home_dir = os.getenv("HOME")
    if not args.stage:
        args.stage = "all"
    if not args.project_name:
        args.project_name = "SAMPLE_"
    if not args.output_dir:
        args.output_dir = os.getcwd()
    if not args.sequences_dir:
        args.sequences_dir = os.getcwd()
    if not args.genomes:
        args.genomes = home_dir + "/ref_genome"
    if not args.genome_version:
        args.genome_version = "hg19"
    if not args.bowtie2:
        args.bowtie2 = "forBowtie2"
    if not args.annotation:
        args.annotation = home_dir + "annotation_file.txt"
    if not args.control_file:
        args.control_file = home_dir + "control_file.txt"
    if not args.refseq_file:
        #TODO: replace this file with newer version!!!
        args.refseq_file = home_dir + "All_genes_data_ct_srt_sed.txt"
    if not args.num_cpus:
        args.num_cpus = "4"
    if not args.ins_annotation:
        args.ins_annotation = home_dir + "ucsc_hg19_ensembl_73_genes_parsed.txt"
    if not args.plot_option:
        args.plot_option = "pdf"
    if not args.fdr_cutoff:
        args.fdr_cutoff = "0.05"

    # set project directory
    project_dir = args.output_dir + "/" + args.project_name

    # create directory structure
    ts.create_output_dir(args.output_dir, args.project_name)
    sub_dir = args.output_dir + "/" + args.project_name
    ts.create_output_dir(sub_dir, "img")
    ts.create_output_dir(sub_dir, "statistics")

    # create log file
    logfile_name = args.output_dir + "/" + args.project_name + "/" \
                   + args.project_name + ".log"
    logging.basicConfig(filename=logfile_name,
                        format='%(levelname)s: %(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.DEBUG)

    # start analysis workflow & logging
    logging.info("Genetic screen workflow 0.0.1")

    # print command line arguments (if debug)
    if args.debug:
        config_param = ts.print_config_param(
            args.project_name,home_dir, args.output_dir,
            args.sequences_dir, project_dir,
            args.sample_file, args.genomes,
            args.genome_version, args.bowtie2,
            str(args.num_cpus)
        )


    # only necessary if --stage is not [all]                   
    sample_file = args.sequences_dir + "/" + args.sample_file

    # load dictionary for file extensions
    data_dir = str(os.path.dirname(os.path.realpath(__file__)).strip('bin')) \
               + "data" + "/"
    print data_dir
    file_ext = ts.load_dictionary(data_dir + 'file_extension.txt')
    # load dictionary for stdout messages
    stdout_msg = ts.load_dictionary(data_dir + 'stdout_message.txt')

    #get execution directory
    dn = os.path.dirname(os.path.realpath(__file__))

    if re.search(r"bam2fastq", args.stage):
        out_filename = project_dir + "/" \
                       + args.project_name + "." \
                       + file_ext['bam2fastq']

        cmd = picard.bam2fastq(
            bamfile=sample_file,
            fastqfile=out_filename
        )

        status = ts.run_cmd(
            message=stdout_msg['bam2fastq'],
            command=cmd,
            debug=args.debug
        )
        sample_file = out_filename

    if re.search(r"all|alignment", args.stage):
        
        out_filename = project_dir + "/" \
                       + args.project_name + "." \
                       + file_ext['alignment']

        cmd = bowtie.alignment(
            genome_path=args.genomes + "/"
                        + args.genome_version,
            fastqfile=sample_file,
            samfile=out_filename,
            num_cpus=str(args.num_cpus),
            metrics=project_dir + "/"
                    + args.project_name
                    + ".align_metrics.txt"
        )

        status = ts.run_cmd(
            message=stdout_msg['alignment'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|filter", args.stage):
        
        sample_file = project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['alignment']
        
        out_filename = project_dir + "/" \
                       + args.project_name + "." \
                       + file_ext['sam2bam']

        cmd = samtools.sam2bam(
            insamfile=sample_file,
            outbamfile=project_dir + "/"
                    + args.project_name + "."
                    + file_ext['sam2bam']
        )

        status = ts.run_cmd(
            message=stdout_msg['sam2bam'],
            command=cmd,
            debug=args.debug
        )

        filter_filename = project_dir + "/" \
                       + args.project_name + "." \
                       + file_ext['filter']

        cmd = samtools.filter_reads(
            inbamfile=out_filename,
            outbamfile=filter_filename,
            mapq="20",
            flag="4"
        )

        status = ts.run_cmd(
            message=stdout_msg['filter'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|duplicates", args.stage):
        
        sample_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['filter']

        cmd = picard.sort_bam(
            inbamfile=sample_file,
            outbamfile=project_dir + "/"
                       + args.project_name + "."
                       + file_ext['sort_bam'],
            sort_order="coordinate"
        )

        status = ts.run_cmd(
            message=stdout_msg['sort_bam'],
            command=cmd,
            debug=args.debug
        )

        sort_bam = project_dir + "/" \
                   + args.project_name + "." \
                   + file_ext['sort_bam']

        cmd = picard.reorder_sam(
            inbamfile=sort_bam,
            outbamfile=project_dir + "/"
                       + args.project_name + "."
                       + file_ext['reorder_sam'],
            genome_path=args.genomes + "/"
                        + args.genome_version
        )

        status = ts.run_cmd(
            message=stdout_msg['reorder_sam'],
            command=cmd,
            debug=args.debug
        )

        reord_bam = project_dir + "/" \
                    + args.project_name + "." \
                    + file_ext['reorder_sam']

        cmd = picard.remove_duplicates(
            inbamfile=reord_bam,
            outbamfile=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['duplicates'],
            metrics_file=project_dir + "/" \
                         + args.project_name + "."
                         + "duplicate_metrics.txt"
        )

        status = ts.run_cmd(
            message=stdout_msg['duplicates'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|index", args.stage):
        
        sample_file = project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['duplicates']

        cmd = samtools.index_bam(
            inbamfile=sample_file
        )

        status = ts.run_cmd(
            message=stdout_msg['index'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|insertions", args.stage):
        
        sample_file = project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['duplicates']
        
        cmd = samtools.bam2sam(
            inbamfile=sample_file,
            outsamfile=project_dir + "/"
                       + args.project_name + "."
                       + file_ext['bam2sam']
        )

        status = ts.run_cmd(
            message=stdout_msg['bam2sam'],
            command=cmd,
            debug=args.debug
        )

        sample_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['bam2sam']

        sam_file_name = project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['insertions']

        sam_out = open(sam_file_name, 'w')
        try:
            sam_out = pf.remove2bpinsertions(
                insamfile=sample_file,
                sam_out=sam_out
            )

        except Exception, exc:
            logging.warning("Error. The reason is: %s", exc)

        finally:
            sam_out.close()

    if re.search(r"all|annotate", args.stage):
        
        sample_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['insertions']

        cmd = samtools.get_header(
            inbamfile=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['sam2bam'],
            header=project_dir + "/"
                       + args.project_name + "."
                       + file_ext['get_header']
        )

        status = ts.run_cmd(
            message=stdout_msg['get_header'],
            command=cmd,
            debug=args.debug
        )

        cmd = ts.concatenate_files(
            file_1=project_dir + "/"
                   + args.project_name + "."
                   + file_ext['get_header'],
            file_2=sample_file,
            output_file=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['cut_header']
        )

        status = ts.run_cmd(
            message=stdout_msg['cut_header'],
            command=cmd,
            debug=args.debug
        )

        sample_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['cut_header']

        cmd = samtools.sam2bam(
            insamfile=sample_file,
            outbamfile=project_dir + "/"
                       + args.project_name + "."
                       + file_ext['header2bam']
        )

        status = ts.run_cmd(
            message=stdout_msg['header2bam'],
            command=cmd,
            debug=args.debug
        )

        sample_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['header2bam']

        cmd = bedtools.sam2bed(
            insamfile=sample_file,
            outbedfile=project_dir + "/"
                        + args.project_name + "."
                        + file_ext['sam2bed']
        )

        status = ts.run_cmd(
            message=stdout_msg['sam2bed'],
            command=cmd,
            debug=args.debug
        )

        sample_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['sam2bed']
        try:

            pf.fix_end_position(
                inbedfile=sample_file,
                outbedfile=project_dir + "/" \
                           + args.project_name + "." \
                           + file_ext['fix_pos']
            )

        except Exception, exc:
            logging.warning("Error. The reason is: %s", exc)

        sample_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['fix_pos']
        
        bed_files = ts.parse_intersectfile(args.annotation)
        # loop through bed files
        for bed in bed_files:
            # annotation_method is only needed for count
            annotation_method = bed[0]
            annotation_name = bed[1]
            annotation_file = bed[2]

            output_file = project_dir + "/" \
                          + args.project_name + "." \
                          + file_ext['fix_pos'] + "." \
                          + annotation_name + ".bed"

            cmd = bedtools.intersectbed(
                inbedfile=sample_file,
                annotation_file=annotation_file,
                outbedfile=output_file
            )

            status = ts.run_cmd(
                message=stdout_msg['fix_pos'],
                command=cmd,
                debug=args.debug
            )

    if re.search(r"all|count", args.stage):

        bed_files = ts.parse_intersectfile(args.annotation)
        # loop through bed files,
        for bed in bed_files:
            # annotation_method is only needed for count
            annotation_method = bed[0]
            annotation_name = bed[1]

            # count insertions
            count_insertions.count_insertions_pythonic(
                inbedfile=project_dir + "/"
                          + args.project_name + "."
                          + file_ext['fix_pos'] + "."
                          + annotation_name + ".bed",
                outfile=project_dir + "/"
                        + args.project_name + "."
                        + annotation_name + "."
                        + file_ext['count'],
                annotation_method=annotation_method,
            )

    if re.search(r"all|fisher", args.stage):

        bed_files = ts.parse_intersectfile(args.annotation)
        # loop through bed files,
        for bed in bed_files:
            annotation_name = bed[1]
            control_file = bed[3]

            cmd = fisher_test.analysis_workflow(
                infile=sample_file,
                control_file=args.control_file,
                outfile=project_dir + "/" \
                        + args.project_name + "." \
                        + file_ext['fisher']
            )

            status = ts.run_cmd(
                message=stdout_msg['fisher'],
                command=cmd,
                debug=args.debug
            )

    if re.search(r"all|plot", args.stage):

        if args.plot_option is "pdf":
            plotWidth = 14
            plotHeight = 10
        else:  # for png
            plotWidth = 3000
            plotHeight = 2400

        sample_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['fisher']

        cmd = rexe.plot_results(
            infile=sample_file,
            refseq_file=args.refseq_file,
            outfile=project_dir + "/" \
                    + args.project_name + "." \
                    + file_ext['bubble'],
            dn=dn
        )

        status = ts.run_cmd(
            message=stdout_msg['bubble'],
            command=cmd,
            debug=args.debug
        )

        insertion_data = project_dir + "/" \
                         + args.project_name + "." \
                         + file_ext['sam2bed']

        cmd = rexe.plot_insertions(
            annotFilePath=args.ins_annotation,
            infile=sample_file,
            insertions=insertion_data,
            outdir=sub_dir + "/" + "img",
            fdrCutoff=args.fdr_cutoff,
            screenName=args.project_name,
            plotOption=args.plot_option,
            plotWidth=plotWidth,
            plotHeight=plotHeight,
            minDistFactor=10,
            dn=dn
        )

        status = ts.run_cmd(
            message=stdout_msg['plot_ins'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|browser", args.stage):

        bed_files = ts.parse_intersectfile(args.annotation)

        for bed in bed_files:

            annotation_name = bed[1]
            
            browser_track.create_track(
                input_file=project_dir + "/" \
                           + args.project_name + "." \
                           + file_ext['fix_pos'] + "." \
                           + annotation_name + ".bed",
                output_file=project_dir + "/" \
                            + args.project_name + "." \
                            + annotation_name + "." \
                            + file_ext['browser'],
                annotation_name=annotation_name
            )

    if re.search(r"all|statistics", args.stage):

        #################################
        print stdout_msg['statistics']  #
        print "Not implemented."        #
        #################################
