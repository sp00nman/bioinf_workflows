#!/usr/bin/env python

# Author: Fiorella Schischlik


import argparse
import re
import logging
from sys import exit
from os import (system, remove, mkdir)
from os.path import (split, splitext, join, exists)
import os
from lib import runnables as rb
from lib import process_files as pf
from utils import tools as ts


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
    parser.add_argument('--project_name', dest='project_name', required=False,
                        help='Name of project directory.')
    parser.add_argument('--output_dir', dest='output_dir', required=False,
                        help='Name of output directory.')
    parser.add_argument('--sequences_dir', dest='sequences_dir', required=False,
                        help='Directory of sequence files.')
    parser.add_argument('--sample_file', dest='sample_file', required=False,
                        help='Input filename.')
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
    parser.add_argument('--refseq_file', dest='refseq_file', required=False,
                        help='Refseq file with start & end position of gene.')
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

    # set project directory
    project_dir = args.output_dir + "/" + args.project_name

    # create directory structure
    ts.create_output_dir(args.output_dir, args.project_name)
    #TODO: create /img and /statistics

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
    data_dir = str(os.path.dirname(os.path.realpath(__file__)).strip('src')) \
               + "data" + "/"
    file_ext = ts.load_dictionary(data_dir + 'file_extension.txt')
    # load dictionary for stdout messages
    stdout_msg = ts.load_dictionary(data_dir + 'stdout_message.txt')

    #get execution directory
    dn = os.path.dirname(os.path.realpath(__file__))

    if re.search(r"bam2fastq", args.stage):
        out_filename = project_dir + "/" \
                       + args.project_name + "." \
                       + file_ext['bam2fastq']

        cmd = rb.bam2fastq(
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

        cmd = rb.alignment(
            genome_path=args.genomes + "/"
                        + args.genome_version,
            fastqfile=sample_file,
            samfile=out_filename,
            num_cpus=str(args.num_cpus),
            metrics=project_dir + "/"
                    + args.project_name
                    + ".align.metrics.txt"
        )

        status = ts.run_cmd(
            message=stdout_msg['bam2fastq'],
            command=cmd,
            debug=args.debug
        )
        sample_file = out_filename

    if re.search(r"all|filter", args.stage):
        out_filename = project_dir + "/" \
                       + args.project_name + "." \
                       + file_ext['sam2bam']

        cmd = rb.sam2bam(
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

        cmd = rb.filter_reads(
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

        sample_file = filter_filename

    if re.search(r"all|duplicates", args.stage):

        cmd = rb.sort_bam(
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

        cmd = rb.reorder_sam(
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

        cmd = rb.remove_duplicates(
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

        sample_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['duplicates']

    if re.search(r"all|index", args.stage):
        
        cmd = rb.index_bam(
            inbamfile=sample_file
        )

        status = ts.run_cmd(
            message=stdout_msg['index'],
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|insertions", args.stage):
        
        cmd = rb.bam2sam(
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
        try:
            sam_out = open(sam_file_name, 'w')

            sam_out = pf.remove2bpinsertions(
                insamfile=sample_file,
                sam_out=sam_out
            )

        except Exception, exc:
            logging.warning("Error. The reason is: %s", exc)

        finally:
            sam_out.close()

        sample_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['insertions']

    if re.search(r"all|annotate", args.stage):

        cmd = rb.get_header(
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

        cmd = rb.concatenate_files(
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

        cmd = rb.sam2bam(
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

        cmd = rb.sam2bed(
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
        
        for bed in bed_files:
            annotation_name = bed[0]
            annotation_file = bed[1]

            output_file = project_dir + "/" \
                          + args.project_name + "." \
                          + file_ext['fix_pos'] + "." \
                          + annotation_name + ".bed"

            cmd = rb.intersectbed(
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

        #TODO: make this more flexible
        exon = project_dir + "/" \
               + args.project_name + "." \
               + file_ext['fix_pos'] + "." \
               + "exon" + "." + "bed"
        intron = project_dir + "/" \
                 + args.project_name + "." \
                 + file_ext['fix_pos'] + "." \
                 + "intron" + "." + "bed"
        output_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['count']

        cmd = pf.count_insertions(
            exonfile=exon,
            intronfile=intron,
            outcountfile=output_file,
            dn=dn
        )

        status = ts.run_cmd(
            message=stdout_msg,
            command=cmd,
            debug=args.debug
        )

        sample_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['count']

    if re.search(r"all|fisher", args.stage):

        output_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['fisher']

        cmd = pf.fisher_test(
            infile=sample_file,
            control_file=args.control_file,
            outfile=project_dir + "/" \
                    + args.project_name + "." \
                    + file_ext['fisher'],
            dn=dn
        )

        status = ts.run_cmd(
            message=stdout_msg,
            command=cmd,
            debug=args.debug
        )

        sample_file = project_dir + "/" \
                      + args.project_name + "." \
                      + file_ext['fisher']

    if re.search(r"all|plot", args.stage):

        cmd = pf.plot_results(
            infile=sample_file,
            refseq_file=args.refseq_file,
            outfile=project_dir + "/" \
                    + args.project_name + "." \
                    + file_ext['bubble'],
            dn=dn
        )

        status = ts.run_cmd(
            message=stdout_msg,
            command=cmd,
            debug=args.debug
        )

    if re.search(r"all|browser", args.stage):

        bed_files = ts.parse_intersectfile(args.annotation)

        for bed in bed_files:
            annotation_name = bed[0]
            pf.browser_track(
                input_file=project_dir + "/" \
                           + args.project_name + "." \
                           + "insertions" + "." \
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
