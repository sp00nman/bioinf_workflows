"""
Collection of functions to execute software for bowtie2
"""


def alignment(genome_path,
              fastqfile,
              samfile,
              num_cpus,
              metrics):
    """
    Alignment of reads to the reference genome with bowtie2.
    :param genome_path: version of the genome, eg [hg19, mm10]
    :param fastqfile: name of FASTQ formatted file
    :param samfile: name of output file
    :param num_cpus: for multi-threading purposes (-p option; bowtie2)
    :return: Command to be executed; type str
    """
    cmd_align = "bowtie2 " \
                "-p %s " \
                "--end-to-end " \
                "--sensitive " \
                "-x %s " \
                "-U %s " \
                "-S %s " \
                "--met-file %s" \
                % (num_cpus,
                   genome_path,
                   fastqfile,
                   samfile,
                   metrics)

    return cmd_align
