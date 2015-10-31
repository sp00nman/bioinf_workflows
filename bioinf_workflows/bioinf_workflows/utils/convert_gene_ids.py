"""
this module maps gene symols or refseq ids to stable ensembl ids
ensgene - ensembl gene id (eg. ENSG....)
genesymb - gene symbols (eg. JAK2)
tab_ensgene_genesymb - hash table containing mapping information
                       for converting gene symbols to ensembl gene ids
tab_ensgene_refseqid - hast table containing
"""

from bioinf_workflows.utils import tools as ts

def read_ensgene_genesymb(filename):
    """
    path = /home/illumina/Fiorella/annotation
    filename = tab_combined_ensgene_genesymb.txt
    #file structure
    ensgene_id      gene_symbol
    ENSG00000000003 TSPAN6
    #key|value key=ensgeneid | value=gene symbol
    """
    ensgene_genesymb = ts.load_tab_delimited(filename)
    tab_ensgene_genesymb = {}

    for line in ensgene_genesymb:
        ensgeneid = line[0]
        genesymb = line[1]
        tab_ensgene_genesymb[genesymb] = ensgeneid
    return tab_ensgene_genesymb


def ensgene2genesymbol(tab_ensgene_genesymb,
                        genesymbol):
    """
    function for converting gene symbols to en
    careful: LOC541471 changed
    """
    ensgeneid = "NOT_FOUND"

    if genesymbol in tab_ensgene_genesymb:
        ensgeneid = tab_ensgene_genesymb[genesymbol]
    else:
        print("ensid not found")
        print genesymbol

    return ensgeneid
