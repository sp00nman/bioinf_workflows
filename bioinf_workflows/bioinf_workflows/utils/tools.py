"""
Collections of utility functions.
"""
import logging
from os import (system, remove, mkdir)
from os.path import (split, splitext, join, exists)


def load_dictionary(file):
    """
    Reads in files of the following format:
    key1: value1
    key2: value2
    :param file: plain text file with
    :return: dictionary of key-value pairs
    """

    d = {}

    try:
        file_handle = open(file)

        for line in file_handle:
            key_value_pair = line.split(':')
            d[key_value_pair[0].strip(' ')] = key_value_pair[1].strip(' ').rstrip('\n')

        file_handle.close()

    except IOError:
        print('Key value read-in file missing.')

    return d


def run_cmd(message, command, debug):
    """
    Print stdout message and log commands, return non-zero exit codes
    if execution failed.
    :param message: Stdout message to be printed.
    :param command: Command to be executed.
    :param debug: Do not execute command if debug is set to 1
    :return: integer, 1 -success, 0 - fail
    """

    print message

    logging.info(message)
    logging.debug(command)
    status = 0

    if not debug:
        status = system(command)
        if status != 0:
            logging.warning("command '%s' returned non-zero "
                            "status: %d'" % (command, status))
    return status


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


def create_output_dir(output_dir,
                      project_name):
    """
    Create project directory.
    """
    if not exists(output_dir + "/" + project_name):
        logging.info('Create folder %s' % output_dir)
        try:
            mkdir(output_dir + "/" + project_name, 0777)
        except IOError, e:
            exit('%s\nFailed to create directory', (e, output_dir))


def load_file(filename):
    """
    Read in files of the following structure
    ID1\tID2\t? --> \t tab separated
    could be 1 or more ids
    :param filename: name of file
    :return: return content of file, type: list
    """

    file_obj = open(filename, 'r')
    try:
        all_content = [line.strip('\n').split('\t') for line in file_obj]
    finally:
        file_obj.close()
    return all_content


def parse_intersectfile(annotation_bed):
    """
    Read in annotation files
    :param annotation_bed: annotation file
    :return: list of BED files
    """

    file_handle = open(annotation_bed)

    try:
        bed_files = [line.rstrip('\r\n').split('\t') \
                     for line in file_handle if line.rstrip('\r\n')]
    except IOError:
        print('Annotation file is missing.')
    file_handle.close()

    return bed_files


def load_tab_delimited(
    input_file,
    sep='\t'
):

    file_handle = open(input_file)

    try:
        vcfs = [line.rstrip('\n').split(sep) for line in file_handle]

    except IOError:
        print('Patient VCF file missing')

    finally:
        file_handle.close()

    return vcfs


def concatenate_files(file_1,
                      file_2,
                      output_file):
    """
    Concatenate two files write output to new file. Uses standard UNIX cat
    :param file_1: 1st parameter
    :param file_2: 2nd parameter
    :param output_file: name of output file
    :return: Command to be executed; type str
    """

    cmd_concatenate_files = "cat %s %s >%s" % (file_1,
                                               file_2,
                                               output_file)
    return cmd_concatenate_files
