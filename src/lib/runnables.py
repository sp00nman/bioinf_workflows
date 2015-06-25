"""
Collection of functions that execute external software.
"""


def bam2fastq(bamfile,
              fastqfile):
    """
    Converts BAM formatted files to FASTQ formatted files with PICARD.
    :param bamfile: Sample filename
    :param fastqfile:project sub-directory name
    :return: Command to be executed; type str
    """
    cmd_bam2fastq = "java -Xmx6g -jar $NGS_PICARD/SamToFastq.jar " \
                    "INPUT=%s " \
                    "FASTQ=%s" % (bamfile, fastqfile)
    return cmd_bam2fastq


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


def sam2bam(samfile,
            bamfile):
    """
    Convert SAM formatted file to BAM formatted file.
    :param samfile: Input samfile
    :param bamfile: Output bamfile
    :return: Command to be executed; type str
    """
    cmd_sam2bam = "samtools view " \
                  "-S " \
                  "-b %s " \
                  ">%s" % (samfile,
                           bamfile)
    return cmd_sam2bam


def filter_reads(inbamfile,
                 outbamfile,
                 mapq,
                 flag):
    """
    Filter read based on
    -q mapping quality of read
    -F only uniquely mapped reads
    :param inbamfile: name of BAM formatted file
    :param outbamfile: name of output file (BAM formatted and filtered)
    :param mapq: mapping quality
    :return: Command to be executed; type str
    """
    cmd_filter = "samtools view " \
                 "-b " \
                 "-q %s " \
                 "-F %s " \
                 "%s >%s" % (mapq,
                             flag,
                             inbamfile,
                             outbamfile)
    return cmd_filter


def sort_bam(project_name,
             project_dir,
             sample_file,
             file_ext):
    """
    Sort sam ? bam file by coordinate.
    :param project_name: name of project (given by user)
    :param output_dir: where the output files should be written
    :return: message to be logged & command to be executed; type str
    """

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "." + file_ext
    msg_sort = "Sort bam file (by coordinate)."
    cmd_sort = "java -Xmx6g -jar $NGS_PICARD/SortSam.jar " \
               "INPUT=%s " \
               "OUTPUT=%s " \
               "SORT_ORDER=coordinate" % (input_file, output_file)
    return msg_sort, cmd_sort


def reorder_sam(project_name,
                project_dir,
                sample_file,
                genomes,
                file_ext):

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "." + file_ext
    msg_reorder = "Reorder bam file."
    cmd_reorder = "java -Xmx6g -jar $NGS_PICARD/ReorderSam.jar " \
                  "INPUT=%s " \
                  "OUTPUT=%s " \
                  "REFERENCE=%s/hg19.fa" % (input_file, output_file, genomes)
    return msg_reorder, cmd_reorder


def count_duplicates(project_name,
                     project_dir,
                     sample_file,
                     file_ext):

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "." + file_ext
    msg_countdup = "Count duplicate reads for each position."
    #TODO: don't hardcode executable files
    cmd_countdup = "~/src/ngsutils/bin/bamutils pcrdup " \
    + "-frag " \
    + "-bam %s " \
    + "-counts %s "
    return msg_countdup, cmd_countdup


def remove_duplicates(project_name,
                      project_dir,
                      sample_file,
                      file_ext):
    """
    Remove duplicate reads.
    :param project_name: name of project (given by user)
    :param output_dir: where the output files should be written
    :return: message to be logged & command to be executed; type str
    """

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "." + file_ext
    metrics_file = project_dir + "/" + project_name + ".duplicates.metrics.txt"
    msg_rmdup = "Remove duplicate reads. "
    cmd_rmdup = "java -Xmx6g -jar $NGS_PICARD/MarkDuplicates.jar " \
                "INPUT=%s " \
                "OUTPUT=%s " \
                "METRICS_FILE=%s " \
                "REMOVE_DUPLICATES=true" % (input_file, output_file,
                                            metrics_file)
    return msg_rmdup, cmd_rmdup


def bam2bai(project_name,
            project_dir,
            sample_file,
            file_ext):

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "." + file_ext
    msg_bam2bai = "Index bam file."
    cmd_bam2bai = "samtools index %s %s" % (input_file, output_file)

    return msg_bam2bai, cmd_bam2bai

def bam2sam(project_name,
            project_dir,
            sample_file,
            file_ext):

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "." + file_ext
    msg_bam2bai = "Convert bam to sam."
    cmd_bam2bai = "samtools view %s > %s" % (input_file, output_file)

    return msg_bam2bai, cmd_bam2bai


def getheader(project_name,
              project_dir,
              file_ext):

    input_file = project_dir + "/" + project_name + ".aligned.bam"
    output_file = project_dir + "/" + project_name + "." + file_ext

    msg_getheader = "Get header."
    cmd_getheader = "samtools view -H %s > %s" % (input_file, output_file)

    return msg_getheader, cmd_getheader


def sam2bed(project_name,
            project_dir,
            sample_file,
            file_ext):

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "." + file_ext
    msg_sam2bed = "Convert samformat to bedformat."
    cmd_sam2bed = "bamToBed -cigar -i %s >%s " % (input_file, output_file)

    return msg_sam2bed, cmd_sam2bed

def intersectbed(project_name,
                 project_dir,
                 sample_file,
                 annotation_file,
                 annotation_name,
                 file_ext):

    input_file = sample_file
    output_file = project_dir + "/" + project_name + "." + file_ext + "." + annotation_name
    msg_intersect = "Intersect " + annotation_file
    cmd_intersect = "intersectBed " \
                    "-a %s " \
                    "-b %s " \
                    "-wo >%s.bed" % (input_file, annotation_file,
                                        output_file)

    return msg_intersect, cmd_intersect