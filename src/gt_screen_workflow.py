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


def config_section_map(section):
    """
    This function was taken from
    https://wiki.python.org/moin/ConfigParserExamples
    :param section:
    :return: dictionary of options
    """
    dict_param = {}
    options = config.options(section)
    for option in options:
        try:
            dict_param[option] = config.get(section, option)
            if dict_param[option] == -1:
                print("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict_param[option] = None
    return dict_param


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
    out_file = project_dir + "/" + project_name + ".filt.aligned.sam"
    msg_filter = "Filter reads with MAPQ< " + mapq + "& non-unique reads."
    cmd_filter = "samtools view -S -q %s -F 4 %s >%s" % (mapq, in_file, out_file)
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
    output_file = project_dir + "/" + project_name + ".sorted.filt.aligned.sam"
    msg_sort = "Sort bam file (by coordinate)."
    cmd_sort = "java -jar $NGS_PICARD/SortSam.jar " \
               "INPUT=%s " \
               "OUTPUT=%s " \
               "SORT_ORDER=coordinate" % (input_file, output_file)
    return msg_sort, cmd_sort


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
    output_file = project_dir + "/" + project_name + ".rm_dupl.sorted.filt.aligned.sam"
    msg_rmdup = "Remove duplicate reads. "
    cmd_rmdup = "java -jar $NGS_PICARD/MarkDuplicates.jar " \
                "INPUT=%s " \
                "OUTPUT=%s " \
                "METRICS_FILE=%s.duplicates.metrics.txt " \
                "REMOVE_DUPLICATES=true" % (input_file, output_file, project_name)

    return msg_rmdup, cmd_rmdup


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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Genetic screen workflow 0.0.1')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False,
                        help='Limit job submission to a particular '
                             'Analysis stage. '
                             '[all,alignment,filter, sort, duplicates,insertions,'
                             'annotate, grouping, count, plot]')
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
    parser.add_argument('--configuration', dest='configuration',
                        required=False, type=str, help='Configuration file (*.ini)')

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

    # set variables
    project_name = args.project_name
    output_dir = args.output_dir
    sequences_dir = args.sequences_dir
    project_dir = output_dir + "/" + project_name
    sample_file = args.sample_file
    genomes = args.genomes
    genome_version = args.genome_version
    bowtie2 = args.bowtie2
    num_cpus = args.num_cpus

    # alternatively read in configuration file
    if args.configuration:
        config = ConfigParser.ConfigParser()
        config.read(args.configuration)

        project_name = config_section_map("project")['project_name']
        home_dir = config_section_map("directories")['home_dir']
        output_dir = config_section_map("directories")['output_dir']
        sequences_dir = config_section_map("directories")['sequences_dir']
        project_dir = config_section_map("directories")['project_dir']
        sample_file = config_section_map("directories")['sample_file']
        genomes = config_section_map("directories")['genomes']
        genome_version = config_section_map("directories")['genome_version']
        bowtie2 = config_section_map("indices")['bowtie2']
        num_cpus = config_section_map("cluster")['num_cpus']


    # create log file
    logfile_name = project_dir + "/" + project_name + ".log"
    logging.basicConfig(filename=logfile_name,
                        format='%(levelname)s: %(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.DEBUG)

    # start analysis workflow & logging
    logging.info("Genetic screen workflow 0.0.1")

    if not args.debug:
        create_output_dir(output_dir, project_name)

    if args.debug:
        config_param = print_config_param(project_name,home_dir,output_dir,sequences_dir,
                           project_dir, sample_file, genomes, genome_version,
                           bowtie2, num_cpus)
        print config_param
        logging.debug(config_param)

    if re.search(r".bam", sample_file):
        (msg, cmd) = bam2fastq(sequences_dir, project_dir,
                               sample_file, project_name, output_dir)
        status = run_cmd(msg, cmd)
        sample_file = project_dir + "/" + project_name + ".fastq"

    if re.search(r"all|alignment", args.stage):
        (msg, cmd) = alignment(genome_version, genomes, project_name,
                               sample_file, project_dir, num_cpus)
        status = run_cmd(msg, cmd)
        sample_file = project_dir + "/" + project_name + ".aligned.sam"

    if re.search(r"all|filter", args.stage):
        (msg, cmd) = filter_reads(project_name, project_dir, sample_file,
                                  mapq="20")
        status = run_cmd(msg, cmd)
        sample_file = project_dir + "/" + project_name + ".filt.aligned.sam"

    if re.search(r"all|duplicates", args.stage):
        (msg, cmd) = sort_bam(project_name, project_dir, sample_file)
        status = run_cmd(msg, cmd)
        sample_file = project_dir + "/" + project_name + ".sorted.filt.aligned.sam"
        (msg, cmd) = remove_duplicates(project_name, output_dir, sample_file)
        status = run_cmd(msg, cmd)



