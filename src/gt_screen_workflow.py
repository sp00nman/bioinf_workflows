#!/usr/bin/env python

# Author: Fiorella Schischlik
# python script to initialize the analysis of genetic screens (eg. genetrap or
# crispr systems in cell lines like KBM7 or BaF3), tested for human and mouse,
# but any other species is feasible

import argparse
import re
import logging
from sys import exit
from os import (system, remove, mkdir)
from os.path import (split, splitext, join, exists)
import pandas as pd
import os
from lib import runnables as rb
from lib.process_files import *
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
                        help='annotation file. Format (tab-separated: exon /path/to/annotation_1.txt \n intron /path/to/annotation_2.txt')
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
        config_param = ts.print_config_param(args.project_name,home_dir, args.output_dir,
                                            args.sequences_dir, project_dir,
                                            args.sample_file, args.genomes, 
                                            args.genome_version, args.bowtie2, 
                                            str(args.num_cpus))


    # only necessary if --stage is not [all]                   
    sample_file = args.sequences_dir + "/" + args.sample_file

    # load dictionary for file extensions
    data_dir = str(os.path.dirname(os.path.realpath(__file__)).strip('src')) + "data" + "/"
    file_ext = ts.load_dictionary(data_dir + 'file_extension.txt')
    # load dictionary for stdout messages
    stdout_msg = ts.load_dictionary(data_dir + 'stdout_message.txt')

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
            samfile=sample_file,
            bamfile=project_dir + "/" \
                    + args.project_name + "." \
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

        (msg, cmd) = sort_bam(args.project_name, project_dir, 
                                sample_file, file_ext=file_ext['sort_bam'])
        status = run_cmd(msg, cmd, args)
        sample_file = project_dir + "/" + args.project_name + "." \
                      + file_ext['sort_bam']
        (msg, cmd) = reorder_sam(args.project_name, project_dir, 
                                    sample_file, args.genomes,
                                    file_ext=file_ext['reorder_sam'])
        status = run_cmd(msg,cmd,args)
        sample_file = project_dir + "/" + args.project_name + "." \
                      + file_ext['reorder_sam']
        (msg, cmd) = remove_duplicates(args.project_name, project_dir,
                                       sample_file, file_ext=file_ext['duplicates'])
        status = run_cmd(msg, cmd,args)
        sample_file = project_dir + "/" + args.project_name + "." \
                      + file_ext['duplicates']

    if re.search(r"all|index", args.stage):
        
        (msg, cmd) = bam2bai(args.project_name, project_dir, 
                                sample_file, file_ext=file_ext['index'])
        status = run_cmd(msg, cmd, args)

    if re.search(r"all|insertions", args.stage):
        
        (msg, cmd) = bam2sam(args.project_name, project_dir, 
                                sample_file, file_ext=file_ext['bam2sam'])
        status = run_cmd(msg, cmd, args)
        sample_file = project_dir + "/" + args.project_name + "." \
                      + file_ext['bam2sam']
        
        #############################
        print stdout_msg['filter']  #
        #############################
        remove2bpinsertions(args.project_name, project_dir, 
                            sample_file, file_ext=file_ext['insertions'])
        sample_file = project_dir + "/" + args.project_name + "." \
                      + file_ext['insertions']

    if re.search(r"all|annotate", args.stage):
        
        #############################
        print stdout_msg['filter']  #
        #############################
        (msg, cmd) = getheader(args.project_name, project_dir,
                               file_ext=file_ext['get_header'])
        status = run_cmd(msg, cmd, args)
        (msg, cmd) = cutheader(args.project_name, project_dir, 
                                sample_file, file_ext=file_ext['cut_header'])
        status = run_cmd(msg, cmd, args)
        sample_file = project_dir + "/" + args.project_name + "." \
                      + file_ext['cut_header']
        (msg, cmd) = sam2bam(args.project_name, project_dir, 
                                sample_file, file_ext=file_ext['header2bam'])
        status = run_cmd(msg, cmd, args)
        sample_file = project_dir + "/" + args.project_name + "." \
                      + file_ext['header2bam']
        (msg, cmd) = sam2bed(args.project_name, project_dir, 
                                sample_file, file_ext=file_ext['sam2bed'])
        status = run_cmd(msg, cmd, args)
        sample_file = project_dir + "/" + args.project_name + "." \
                      + file_ext['sam2bed']
        
        ##############################
        print stdout_msg['fix_pos']  #
        ##############################
        fix_end_position(args.project_name, project_dir,
                            sample_file, file_ext=file_ext['fix_pos'])
        
        sample_file = project_dir + "/" + args.project_name + "." \
                      + file_ext['fix_pos']
        
        #################################
        print stdout_msg['insertions']  #
        #################################
        
        bed_files = parse_intersectfile(args.annotation)
        
        for bed in bed_files:
            annotation_name = bed[0]
            annotation_file = bed[1]
        
            (msg, cmd) = intersectbed(args.project_name, project_dir,
                                    sample_file, annotation_file, 
                                    annotation_name=annotation_name,
                                    file_ext="insertions" )
            status = run_cmd(msg, cmd, args)

    if re.search(r"all|count", args.stage):
        
        ############################
        print stdout_msg['count']  #
        ############################
        (msg, cmd) = count_insertions(args.project_name, project_dir,
                                      file_ext=file_ext['count'])
        status = run_cmd(msg, cmd, args)
        sample_file = project_dir + "/" + args.project_name + "." \
                      + file_ext['count']

    if re.search(r"all|fisher", args.stage):
        
        #############################
        print stdout_msg['filter']  #
        #############################
        (msg, cmd) = fisher_test(args.project_name, project_dir, sample_file,
                                 args.control_file, file_ext=file_ext['fisher'])
        status = run_cmd(msg, cmd, args)
        sample_file = project_dir + "/" + args.project_name + "." \
                      + file_ext['fisher']

    if re.search(r"all|plot", args.stage):
        
        #############################
        print stdout_msg['bubble']  #
        #############################
        (msg, cmd) = plot_results(args.project_name, project_dir,
                                  args.refseq_file, sample_file,
                                  file_ext=file_ext['bubble'])
        status = run_cmd(msg, cmd, args)

    if re.search(r"all|browser", args.stage):

        ##############################
        print stdout_msg['browser']  #
        ##############################
        bed_files = parse_intersectfile(args.annotation)

        for bed in bed_files:
            annotation_name = bed[0]
            browser_track(args.project_name, project_dir,
                          annotation_name=annotation_name,
                          file_ext=file_ext['browser'])

    if re.search(r"all|statistics", args.stage):

        #################################
        print stdout_msg['statistics']  #
        print "Not implemented."        #
        #################################





