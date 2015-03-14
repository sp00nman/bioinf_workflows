#!/usr/bin/env python

# Author: Fiorella Schischlik
# python script to initialize the analysis of genetic screens (eg. genetrap or
# crispr systems in cell lines like KBM7 or BaF3), tested for human and mouse,
# but any other species is feasible

import argparse
import re
import ConfigParser
import logging
from sys import exit
from os import (system, remove, mkdir)
from os.path import (split, splitext, join, exists)
import os
import pysam

def print_config_param(project_name,
                       home_dir,
                       output_dir,
                       sequences_dir,
                       project_dir,
                       sample_file,
                       genomes,
                       genome_version,
                       bowtie2,
                       num_cpus):

    config_param = "[project name:" + project_name + ", " \
                   + "home directory:" + home_dir + ", " \
                   + "output directory:" + output_dir + ", " \
                   + "sequence directory:" + sequences_dir + ", " \
                   + "project directory:" + project_dir + ", " \
                   + "sample file name:" + sample_file + ", " \
                   + "genome directory:" + genomes + ", " \
                   + "genome version:" + genome_version + ", " \
                   + "indices of bowtie2:" + bowtie2 + ", " \
                   + "number of cpus:" + num_cpus + "]"

    return config_param


def run_cmd(msg, cmd):
    logging.info(msg)
    logging.debug(cmd)
    status = 0
    if not args.debug:
        status = system(cmd)
        if status != 0:
            logging.warning("command '%s' returned non-zero "
                            "status: %d'" % (cmd, status))
    return status


def create_output_dir(output_dir,
                      project_name):
    """
    :param output_dir:
    :param project_name:
    :return:
    """
    if not exists(output_dir + "/" + project_name):
        logging.info('Create folder %s' % output_dir)
        try:
            mkdir(output_dir + "/" + project_name, 0777)
        except IOError, e:
            exit('%s\nFailed to create directory', (e, output_dir))


def bam2fastq(sequences_dir,
              project_dir,
              sample_file,
              project_name,
              output_dir):
    """
    Converts bam file to fastq files with PICARD
    :param sequences_dir: Sub-directory for processed sequences in
                          Illumina2bam format
    :param project_dir:project sub-directory name
    :param sample_file: sample file name
    :param project_name: name of the project (given by user)
    :return: message to be logged & command to be executed; type str
    """

    input_file = sequences_dir + "/" + sample_file
    msg_bam2fastq = "Convert bam to fastq."
    cmd_bam2fastq = "java -jar $NGS_PICARD/SamToFastq.jar " \
                    "INPUT=%s " \
                    "FASTQ=%s/%s.fastq" % (input_file, project_dir, project_name)
    return msg_bam2fastq, cmd_bam2fastq


def alignment(genome_version,
              genomes,
              project_name,
              sample_file,
              project_dir,
              num_cpus):

    """
    Aligns fastq reads to the reference genome with bowtie2
    :param genome_version: version of the genome, eg [hg19, mm10]
    :param genomes: path of the genomes file
    :param project_name: name of the project (given by user)
    :param sequences_dir: Sub-directory for processed sequences in
                          Illumina2bam format
    :param project_dir:project sub-directory name
    :param sample_file: sample file name
    :param output_dir: where the output files should be written
    :param num_cpus: for multi-threading purposes (-p option; bowtie2)
    :return: message to be logged & command to be executed; type str
    """

    genome_path = genomes + "/" + genome_version
    sample_path_file = sample_file
    out_file = project_dir + "/" + project_name + ".aligned.sam"
    out_file_metrics = project_dir + "/" + project_name + ".align.metrics.txt"

    msg_align = "Mapping reads to genome " + genome_version
    cmd_align = "bowtie2 -p %s --end-to-end --sensitive -x %s -U %s " \
                "-S %s --met-file %s" \
                % (num_cpus, genome_path, sample_path_file,
                   out_file, out_file_metrics)
    return msg_align, cmd_align


def sam2bam(project_name,
            project_dir,
            sample_file):

    in_file = sample_file
    out_file = project_dir + "/" + project_name + ".bam"
    msg_sam2bam = "Convert sam to bam format."
    cmd_sam2bam = "samtools view -S -b %s >%s" % (in_file, out_file)
    return msg_sam2bam, cmd_sam2bam


def filter_reads(project_name,
                 project_dir,
                 sample_file,
                 mapq):
    """
    Filter read based on
    :param project_name: name of the project (given by user)
    :param output_dir: where the output files should be written
    :param mapq: mapping quality, user defined
    :return: message to be logged & command to be executed; type str
    """

    in_file = sample_file
    out_file = project_dir + "/" + project_name + ".filt.aligned.bam"
    msg_filter = "Filter reads with MAPQ< " + mapq + "& non-unique reads."
    cmd_filter = "samtools view -b -q %s -F 4 %s >%s" % (mapq, in_file, out_file)
    return msg_filter, cmd_filter


def sort_bam(project_name,
             project_dir,
             sample_file):
    """
    Sort sam ? bam file by coordinate.
    :param project_name: name of project (given by user)
    :param output_dir: where the output files should be written
    :return: message to be logged & command to be executed; type str
    """

    input_file = sample_file
    output_file = project_dir + "/" + project_name + ".sorted.filt.aligned.bam"
    msg_sort = "Sort bam file (by coordinate)."
    cmd_sort = "java -jar $NGS_PICARD/SortSam.jar " \
               "INPUT=%s " \
               "OUTPUT=%s " \
               "SORT_ORDER=coordinate" % (input_file, output_file)
    return msg_sort, cmd_sort


def reorder_sam(project_name,
                project_dir,
                sample_file,
                genomes):

    input_file = sample_file
    output_file = project_dir + "/" + project_name + ".reorder.bam"
    msg_reorder = "Reorder bam file."
    cmd_reorder = "java -jar $NGS_PICARD/ReorderSam.jar " \
                  "INPUT=%s " \
                  "OUTPUT=%s " \
                  "REFERENCE=%s/hg19.fa" % (input_file, output_file, genomes)
    return msg_reorder, cmd_reorder


def count_duplicates(project_name,
                     project_dir,
                     sample_file):

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "."
    msg_countdup = "Count duplicate reads for each position."
    #TODO: don't hardcode executable files
    cmd_countdup = "~/src/ngsutils/bin/bamutils pcrdup " \
    "-frag " \
    "-bam %s " \
    "-counts %s " \


def remove_duplicates(project_name,
                      project_dir,
                      sample_file):
    """
    Remove duplicate reads.
    :param project_name: name of project (given by user)
    :param output_dir: where the output files should be written
    :return: message to be logged & command to be executed; type str
    """

    input_file = sample_file
    output_file = project_dir + "/" + project_name + ".rm_dupl.sorted.filt.aligned.bam"
    metrics_file = project_dir + "/" + project_name + ".duplicates.metrics.txt"
    msg_rmdup = "Remove duplicate reads. "
    cmd_rmdup = "java -jar $NGS_PICARD/MarkDuplicates.jar " \
                "INPUT=%s " \
                "OUTPUT=%s " \
                "METRICS_FILE=%s " \
                "REMOVE_DUPLICATES=true" % (input_file, output_file,
                                            metrics_file)
    return msg_rmdup, cmd_rmdup


def bam2bai(project_name,
            project_dir,
            sample_file):

    input_file = sample_file
    output_file = project_dir + "/" + project_name \
                  + ".rm_dupl.sorted.filt.aligned.bam.bai"
    msg_bam2bai = "Index bam file."
    cmd_bam2bai = "samtools index %s %s" % (input_file, output_file)

    return msg_bam2bai, cmd_bam2bai

def bam2sam(project_name,
            project_dir,
            sample_file):

    input_file = sample_file
    output_file = project_dir + "/" + project_name \
                  + ".rm_dupl.sorted.filt.aligned.sam"
    msg_bam2bai = "Convert bam to sam."
    cmd_bam2bai = "samtools view %s > %s" % (input_file, output_file)

    return msg_bam2bai, cmd_bam2bai


def load_files(filename):
    """
    read in files of the following structure
    ID1\tID2\t? --> \t tab separated
    could be 1 or more ids
    """
    file_obj = open(filename, 'r')
    try:
        all_content = [line.strip('\n').split('\t') for line in file_obj]
    finally:
        file_obj.close()
    return all_content


def remove2bpinsertions(project_name,
                        project_dir,
                        sample_file):

    input_file = sample_file
    output_file = project_dir + "/" + project_name \
                  + ".rm2bp_insertions.sam"

    sam_file = load_files(input_file)
    sam_out = open(output_file, 'w')

    for index in range(len(sam_file)):
        print index
        if index == 0:
            print '\t'.join(sam_file[index])
            sam_out.writelines('\t'.join(sam_file[index-1]))
            sam_out.writelines("\n")
            continue

        # samfile chromosome position is column 4
        position = int(sam_file[index-1][3])

        next_read = sam_file[index]
        next_position = int(next_read[3])

        # bpdis --> base pair distance
        bp_dist = next_position-position

        if not (bp_dist <= 2):
            print '\t'.join(sam_file[index])
            sam_out.writelines('\t'.join(sam_file[index]))
            sam_out.writelines("\n")

    sam_out.close()
    return None


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Genetic screen workflow 0.0.1')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False,
                        help='Limit job submission to a particular '
                             'Analysis stage. '
                             '[all,alignment,filter, sort, duplicates, index,'
                             'insertions, annotate, grouping, count, plot]')
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
    if not args.num_cpus:
        args.num_cpus = "4"

    # set project directory
    project_dir = args.output_dir + "/" + args.project_name

    # creat project directory
    create_output_dir(args.output_dir, args.project_name)

    # create log file
    logfile_name = project_dir + "/" + args.project_name + ".log"
    logging.basicConfig(filename=logfile_name,
                        format='%(levelname)s: %(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.DEBUG)

    # start analysis workflow & logging
    logging.info("Genetic screen workflow 0.0.1")

    # print command line arguments (if debug)
    if args.debug:
        config_param = print_config_param(args.project_name,home_dir, args.output_dir,
                                            args.sequences_dir, project_dir, 
                                            args.sample_file, args.genomes, 
                                            args.genome_version, args.bowtie2, 
                                            str(args.num_cpus)) 
    # start workflow
    if re.search(r".bam", args.sample_file):
        (msg, cmd) = bam2fastq(args.sequences_dir, project_dir,
                               args.sample_file, args.project_name, 
                               args.output_dir)
        status = run_cmd(msg, cmd)

    if re.search(r"all|alignment", args.stage):
        if re.search(r".fastq", args.sample_file):
            sample_file = args.sequences_dir + "/" + args.sample_file
        else:
            sample_file = project_dir + "/" + args.project_name + ".fastq"

        (msg, cmd) = alignment(args.genome_version, args.genomes, 
                                args.project_name, sample_file, 
                                project_dir, str(args.num_cpus))
        status = run_cmd(msg, cmd)
        sample_file = project_dir + "/" + args.project_name + ".aligned.sam"

    if re.search(r"all|filter", args.stage):
        (msg, cmd) = sam2bam(args.project_name, project_dir, sample_file)
        status = run_cmd(msg,cmd)
        sample_file = project_dir + "/" + args.project_name + ".bam"
        (msg, cmd) = filter_reads(args.project_name, project_dir, sample_file,
                                  mapq="20")
        status = run_cmd(msg, cmd)
        sample_file = project_dir + "/" + args.project_name + ".filt.aligned.bam"

    if re.search(r"all|duplicates", args.stage):
        (msg, cmd) = sort_bam(args.project_name, project_dir, sample_file)
        status = run_cmd(msg, cmd)
        sample_file = project_dir + "/" + args.project_name + ".sorted.filt.aligned.bam"
        (msg, cmd) = reorder_sam(args.project_name, project_dir, sample_file, args.genomes)
        status = run_cmd(msg,cmd)
        sample_file = project_dir + "/" + args.project_name + ".reorder.bam"
        (msg, cmd) = remove_duplicates(args.project_name, project_dir, sample_file)
        status = run_cmd(msg, cmd)

    if re.search(r"all|index", args.stage):
        (msg, cmd) = bam2bai(args.project_name, project_dir, sample_file)
        status = run_cmd(msg, cmd)
        sample_file = project_dir + "/" + args.project_name + ".sorted.filt.aligned.bam"

    if re.search(r"all|bam2sam", args.stage):
        (msg, cmd) = bam2sam(args.project_name, project_dir, sample_file)
        status = run_cmd(msg, cmd)
        sample_file = project_dir + "/" + args.project_name + ".rm_dupl.sorted.filt.aligned.sam"

    if re.search(r"all|insertions", args.stage):
        (msg, cmd) = remove2bpinsertions(args.project_name, project_dir, sample_file)
        status = run_cmd(msg, cmd)
        sample_file = project_dir + "/" + args.project_name + ".rm2bp_insertions.sam"
