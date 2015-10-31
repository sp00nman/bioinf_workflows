"""
Collection of functions to execute R scripts located in bioinf_workflows/bin
"""


def count_insertions_r(exonfile,
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


def fisher_test_r(infile,
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


def plot_insertions(annotFilePath,
                    infile,
                    insertions,
                    outdir,
                    fdrCutoff,
                    screenName,
                    plotOption,
                    plotWidth,
                    plotHeight,
                    minDistFactor,
                    dn):
    """
    Plot insertions for gene above a specified FDR cutoff.
    """
    cmd_plot_insertions = "Rscript --vanilla " + dn \
                          + "/" + "create_insertion_plots_150721.R " \
                          + "%s %s %s %s %s %s %s %s %s %s " % \
                            (annotFilePath, infile, insertions, outdir,
                             fdrCutoff, screenName, plotOption,
                             plotWidth,plotHeight, minDistFactor)

    return cmd_plot_insertions




