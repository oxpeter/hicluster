#!/usr/bin/env python

""" collects gene lists and works out all common pathways/GO terms for sets of lists

"""
import string
import time
import sys, os, re
import argparse

import numpy
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import scipy
from scipy import stats

from genomepy import genematch, config
from matplotlib_venn import venn2, venn3, venn2_unweighted, venn3_unweighted
import hicluster

########################################################################################

def define_arguments():
    parser = argparse.ArgumentParser(description="Performs set analysis on multiple lists and outputs venn diagrams.")

    # input options
    parser.add_argument("experiments", metavar='filename', nargs='+', type=str, help="a data file for comparing")
    parser.add_argument("-B", "--background", type=str, help="Provide a file for using as background")
    parser.add_argument("-G", "--genes_only", action='store_true',  help="Turn off pathway comparison (much much faster)")
    parser.add_argument("-v", "--verbose", action='store_true',default=False,  help="print messages")
    parser.add_argument("-S", "--show_common", action='store_true',  help="Show all common genes")
    parser.add_argument("-U", "--unweighted", action='store_true',  help="draw unweighted circles")
    parser.add_argument("-C", "--column", type=int, default=0, help="the column to extract data from. (default = 0)")

    return parser

def run_arguments(args):

    verbalise = hicluster.check_verbose(args.verbose)


    experiments = {} # dictionary to contain the lists of genes
    i = 0
    for file in args.experiments:
        experiments[str(i) + '-' + os.path.basename(file)] = hicluster.make_a_list(file, args.column)
        i += 1

    # check correct number of files supplied:
    if  len(experiments) > 4 or len(experiments) < 2:
        verbalise("R", "can only compare two to four experiments. You supplied %d!" % (len(experiments)))
        exit()

    go_obj = genematch.GO_maker()


    gene_pathways = {}
    gene_sets = {}
    path_sets = {}

    for exp in experiments:
        if args.genes_only:
            gene_sets['G' + exp] = set(experiments[exp])
        else:
            verbalise("Y", "Collecting paths for geneset", exp )

            # collect pathways and GO terms for each gene:
            gene_pathways[exp]= genematch.cbir_to_pathway(experiments[exp]).values()

            for godict in go_obj.fetch_gos(experiments[exp]):
                gene_pathways[exp] += godict.keys()

            # error checking to make sure only strings have been sent to the dictionary:
            for pathway in gene_pathways[exp]:
                if isinstance(pathway, list):
                    print pathway

            # convert gene and pathway lists to sets:
            gene_sets['G' + exp] = set(experiments[exp])
            path_sets['P' + exp] = set(gene_pathways[exp])


    if args.genes_only:
        if len(experiments) == 3:
            venn_3way([gene_sets[exp] for exp in gene_sets], [gene_sets[exp] for exp in gene_sets], [exp for exp in gene_sets], show_common=args.show_common, unweighted=args.unweighted)
        elif len(experiments) == 2:
            venn_2way([gene_sets[exp] for exp in gene_sets], [gene_sets[exp] for exp in gene_sets], [exp for exp in gene_sets], show_common=args.show_common, unweighted=args.unweighted)
            set1, set2 = [gene_sets[exp] for exp in gene_sets]

            # GO enrichment of common genes:
            pvals = genematch.go_enrichment(set1 & set2)
            if len(pvals) > 2:
                qvals = hicluster.p_to_q(pvals.values(), display_on=True, cut1s=True)
                #print [qvals[p] for p in qvals if qvals[p] < 0.1]

                background = hicluster.make_a_list(args.background)
                screened = [gene for gene in background if gene not in set1 & set2]

                # kegg enrichment of common genes
                pvalupper, pvalpath = genematch.kegg_pathway_enrichment(set1 & set2, screened)

                qvalspath = hicluster.p_to_q(pvalpath.values(), display_on=True, cut1s=False, conservative=True)
                qvalsupper = hicluster.p_to_q(pvalupper.values(), display_on=True, cut1s=False, conservative=True)

                verbalise( "G","### q values of pathways ###")
                verbalise( "\n".join([p + " " + str(qvalspath[p]) for p in qvalspath if qvalspath[p] < 0.1]))
                verbalise("G", "### q values of pathway types ###")
                verbalise( "\n".join([p + " " + str(qvalsupper[p]) for p in qvalsupper if qvalsupper[p] < 0.1]))


                verbalise( "G","### p values of pathways ###")
                verbalise( "\n".join([p + " " + str(pvalpath[p]) for p in pvalpath if pvalpath[p] < 0.05]))
                verbalise( "G","### p values of pathway types ###")
                verbalise( "\n".join([p + " " + str(pvalupper[p]) for p in pvalupper if pvalupper[p] < 0.05]))

        elif len(experiments) == 4:
            venn_4by4([gene_sets[exp] for exp in gene_sets], [gene_sets[exp] for exp in gene_sets], [exp for exp in gene_sets], show_common=args.show_common, unweighted=args.unweighted)

    else:
        if len(experiments) == 3:
            venn_3way([gene_sets[exp] for exp in gene_sets], [path_sets[exp] for exp in path_sets], [exp for exp in gene_sets], show_common=args.show_common, unweighted=args.unweighted)
        elif len(experiments) == 2:
            venn_2way([gene_sets[exp] for exp in gene_sets], [path_sets[exp] for exp in path_sets], [exp for exp in gene_sets], show_common=args.show_common, unweighted=args.unweighted)
        elif len(experiments) == 4:
            venn_4by4([gene_sets[exp] for exp in gene_sets], [path_sets[exp] for exp in path_sets], [exp for exp in gene_sets], show_common=args.show_common, unweighted=args.unweighted)

def venn_3way(genesets, pathsets, names, show_common=False, unweighted=False):
    if unweighted:
        draw = venn3_unweighted
    else:
        draw = venn3

    figure, axes = plt.subplots(1, 2)
    v1 = draw([genesets[0], genesets[1], genesets[2]], (names[0], names[1], names[2]), ax=axes[0])
    v2 = draw([pathsets[0], pathsets[1], pathsets][2], (names[0], names[1], names[2]), ax=axes[1])
    plt.show()
    if show_common:
        for gene in genesets[0] & genesets[1] & genesets[2]:
            print gene

def venn_4by4(genesets, pathsets, names, show_common=False, unweighted=False):

    if unweighted:
        draw = venn2_unweighted
    else:
        draw = venn2

    figure, axes = plt.subplots(2, 3)

    figure.text(0.45, 0.95, 'Overlap of Orthologous Differentially Expressed Genes', ha="center", va="bottom", size="large")
    figure.text(0.05, 0.9, names[0], ha="left", va="bottom", size="medium",color="yellow")
    figure.text(0.05, 0.875, names[1], ha="left", va="bottom", size="medium", color="cyan")
    figure.text(0.05,0.85,names[2], ha="left", va="bottom", size="medium",color="green")
    figure.text(0.05,0.825,names[3], ha="left", va="bottom", size="medium",color="red")
    v1 = draw([genesets[0], genesets[1]], set_labels = (names[0][:2],names[1][:2]), ax=axes[0][0])
    v2 = draw([genesets[0], genesets[2]], set_labels = (names[0][:2],names[2][:2]), ax=axes[0][1])
    v3 = draw([genesets[0], genesets[3]], set_labels = (names[0][:2],names[3][:2]), ax=axes[0][2])
    v4 = draw([genesets[1], genesets[2]], set_labels = (names[1][:2],names[2][:2]), ax=axes[1][0])
    v5 = draw([genesets[1], genesets[3]], set_labels = (names[1][:2],names[3][:2]), ax=axes[1][1])
    v6 = draw([genesets[2], genesets[3]], set_labels = (names[2][:2],names[3][:2]), ax=axes[1][2])
    v1.get_patch_by_id('10').set_color('yellow')
    v2.get_patch_by_id('10').set_color('yellow')
    v3.get_patch_by_id('10').set_color('yellow')
    v1.get_patch_by_id('01').set_color('cyan')
    v2.get_patch_by_id('01').set_color('green')
    v3.get_patch_by_id('01').set_color('red')

    v4.get_patch_by_id('10').set_color('cyan')
    v5.get_patch_by_id('10').set_color('cyan')
    v6.get_patch_by_id('10').set_color('green')
    v4.get_patch_by_id('01').set_color('green')
    v5.get_patch_by_id('01').set_color('red')
    v6.get_patch_by_id('01').set_color('red')
    try:
        v4.get_patch_by_id('11').set_color('brown')
        v5.get_patch_by_id('11').set_color('blue')
        v6.get_patch_by_id('11').set_color('brown')
        v1.get_patch_by_id('11').set_color('green')
        v2.get_patch_by_id('11').set_color('brown')
        v3.get_patch_by_id('11').set_color('orange')
    except AttributeError:
        pass

    #plt.annotate('Unknown effect', xy=v2.get_label_by_id('100').get_position() - numpy.array([0, 0.05]), xytext=(-70,-70),
    #        ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    #        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    plt.show()

    figure, axes = plt.subplots(2, 3)

    figure.text(0.45, 0.95, "Overlap of Orthologous Differentially Expressed Genes' pathways", ha="center", va="bottom", size="large")
    figure.text(0.05, 0.9, names[0], ha="left", va="bottom", size="medium",color="yellow")
    figure.text(0.05, 0.875, names[1], ha="left", va="bottom", size="medium", color="cyan")
    figure.text(0.05,0.85,names[2], ha="left", va="bottom", size="medium",color="green")
    figure.text(0.05,0.825,names[3], ha="left", va="bottom", size="medium",color="red")

    v1 = draw([pathsets[0], pathsets[1]], set_labels = (names[0][:2],names[1][:2]), ax=axes[0][0])
    v2 = draw([pathsets[0], pathsets[2]], set_labels = (names[0][:2],names[2][:2]), ax=axes[0][1])
    v3 = draw([pathsets[0], pathsets[3]], set_labels = (names[0][:2],names[3][:2]), ax=axes[0][2])
    v4 = draw([pathsets[1], pathsets[2]], set_labels = (names[1][:2],names[2][:2]), ax=axes[1][0])
    v5 = draw([pathsets[1], pathsets[3]], set_labels = (names[1][:2],names[3][:2]), ax=axes[1][1])
    v6 = draw([pathsets[2], pathsets[3]], set_labels = (names[2][:2],names[3][:2]), ax=axes[1][2])
    v1.get_patch_by_id('10').set_color('yellow')
    v2.get_patch_by_id('10').set_color('yellow')
    v3.get_patch_by_id('10').set_color('yellow')
    v1.get_patch_by_id('01').set_color('cyan')
    v2.get_patch_by_id('01').set_color('green')
    v3.get_patch_by_id('01').set_color('red')

    v4.get_patch_by_id('10').set_color('cyan')
    v5.get_patch_by_id('10').set_color('cyan')
    v6.get_patch_by_id('10').set_color('green')
    v4.get_patch_by_id('01').set_color('green')
    v5.get_patch_by_id('01').set_color('red')
    v6.get_patch_by_id('01').set_color('red')
    try:
        v4.get_patch_by_id('11').set_color('brown')
        v5.get_patch_by_id('11').set_color('blue')
        v6.get_patch_by_id('11').set_color('brown')
        v1.get_patch_by_id('11').set_color('green')
        v2.get_patch_by_id('11').set_color('brown')
        v3.get_patch_by_id('11').set_color('orange')
    except AttributeError:
        pass

    plt.show()

def venn_2way(genesets, pathsets, names, show_common=False, unweighted=False):
    if unweighted:
        draw = venn2_unweighted
    else:
        draw = venn2

    figure, axes = plt.subplots(1,2)
    figure.text(0.45, 0.95, "Overlap of Orthologous Differentially Expressed Genes", ha="center", va="bottom", size="large")
    figure.text(0.05, 0.9, names[0], ha="left", va="bottom", size="medium",color="red")
    figure.text(0.05, 0.875, names[1], ha="left", va="bottom", size="medium", color="cyan")

    v1 = draw([genesets[0], genesets[1]], set_labels = (names[0][:2],names[1][:2]), ax=axes[0])
    v2 = draw([pathsets[0], pathsets[1]], set_labels = (names[0][:2],names[1][:2]), ax=axes[1])

    v1.get_patch_by_id('10').set_color('red')
    v2.get_patch_by_id('10').set_color('red')
    v1.get_patch_by_id('01').set_color('cyan')
    v2.get_patch_by_id('01').set_color('cyan')
    try:
        v1.get_patch_by_id('11').set_color('blue')
        v2.get_patch_by_id('11').set_color('blue')
    except AttributeError:
        pass
    plt.show()

    if show_common:
        dbpaths = config.import_paths()
        ncbi = genematch.cbir_ncbi(list(genesets[0] & genesets[1]), dbpaths=dbpaths)

        for gene in ncbi :
            if ncbi[gene] == "no NCBI item":
                print gene
            else:
                print ncbi[gene]



def old():
    figure, axes = plt.subplots(2, 2)

    venn2(subsets=(1090,311,39), set_labels = ('Cerapachys', 'Polistes'), ax=axes[0][0])
    venn2(subsets=(1030,827,99), set_labels = ('Cerapachys', 'Repro W / Q'), ax=axes[0][1])
    venn2(subsets=(1025,850,104), set_labels = ('Cerapachys', 'Sterile W / Q'), ax=axes[1][0])
    venn2(subsets=(1122,76,7), set_labels = ('Cerapachys', 'Repro W / Sterile W'), ax=axes[1][1])
    plt.show()

    figure, axes = plt.subplots(1,2)
    venn3(subsets = (5, 3, 5, 2, 1, 2, 2), set_labels=('Polistes', 'Cerapachys', 'Apis'), ax=axes[0])
    venn3(subsets = (44, 35, 46, 12, 11, 2, 10), set_labels=('Polistes', 'Cerapachys', 'Apis'), ax=axes[1])
    plt.show()


    experiments = {}
    experiments['RW-Q'] = '/Volumes/antqueen/genomics/experiments/analyses/PRO20140917_broodswap_controls/comparative/Gro.repro_work-queen.comp.list'
    experiments['SW-Q'] = '/Volumes/antqueen/genomics/experiments/analyses/PRO20140917_broodswap_controls/comparative/Gro.sterile_work-queen.comp.list'
    experiments['SW-RW'] = '/Volumes/antqueen/genomics/experiments/analyses/PRO20140917_broodswap_controls/comparative/Gro.sterile_work-repro_worker.comp.list'
    experiments['Polistes'] = '/Volumes/antqueen/genomics/experiments/analyses/PRO20140917_broodswap_controls/comparative/Toth.Polistes.comp.list'
    experiments['Pol_QW'] = '/Volumes/antqueen/genomics/experiments/analyses/PRO20140917_broodswap_controls/comparative/P.metricus_QvsW.Cbir_orthologs.txt'
    experiments['cbir_meth'] = ''
    experiments['cbir_bs'] = '/Volumes/antqueen/genomics/experiments/analyses/PRO20140917_broodswap_controls/clustering/clustering.Sep22_18.28.t_test.list'

    print "Making GO class"
    go_obj = genematch.GO_maker()

    diflists = {}
    gene_pathways = {}
    pathway_sets = {}

    for exp in experiments:
        print "Collecting paths for geneset", exp
        diflists[exp] = hicluster.make_a_list(experiments[exp])

        pathway_sets['G' + exp] = set(diflists[exp])

        gene_pathways[exp]= genematch.cbir_to_pathway(diflists[exp]).values()

        for godict in go_obj.fetch_gos(diflists[exp]):
            gene_pathways[exp] += godict.keys()
        for pathway in gene_pathways[exp]:
            if isinstance(pathway, list):
                print pathway
        pathway_sets[exp] = set(gene_pathways[exp])

    # should now have a dictionary with a set of kegg and GO pathways
    # for each RNA_brain experiment.
    pathway_sets['GallCbir']   = pathway_sets['Gcbir_meth'] & pathway_sets['Gcbir_bs']
    pathway_sets['Gpol-cbir']  = pathway_sets['GPolistes']  & pathway_sets['Gcbir_bs']
    pathway_sets['GRWQ-cbir']  = pathway_sets['GRW-Q']  & pathway_sets['Gcbir_bs']
    pathway_sets['GSWQ-cbir']  = pathway_sets['GSW-Q']  & pathway_sets['Gcbir_bs']
    pathway_sets['GSWRW-cbir'] = pathway_sets['GSW-RW'] & pathway_sets['Gcbir_bs']
    pathway_sets['Gpol_QW-cbir'] = pathway_sets['GPol_QW'] & pathway_sets['Gcbir_bs']

    print "gene sets:"
    print "(Cbir_BS %d (%d) %d Cbir_meth)" % (len(pathway_sets['Gcbir_bs']) - len(pathway_sets['GallCbir']), len(pathway_sets['GallCbir']), len(pathway_sets['Gcbir_meth']) - len(pathway_sets['GallCbir']))
    print "(Cbir_BS %d (%d) %d Polistes)" % (len(pathway_sets['Gcbir_bs']) - len(pathway_sets['Gpol-cbir']), len(pathway_sets['Gpol-cbir']), len(pathway_sets['GPolistes']) - len(pathway_sets['Gpol-cbir']))
    print "(Cbir_BS %d (%d) %d RW-Q)" % (len(pathway_sets['Gcbir_bs']) - len(pathway_sets['GRWQ-cbir']), len(pathway_sets['GRWQ-cbir']), len(pathway_sets['GRW-Q']) - len(pathway_sets['GRWQ-cbir']))
    print "(Cbir_BS %d (%d) %d SW-Q)" % (len(pathway_sets['Gcbir_bs']) - len(pathway_sets['GSWQ-cbir']), len(pathway_sets['GSWQ-cbir']), len(pathway_sets['GSW-Q']) - len(pathway_sets['GSWQ-cbir']))
    print "(Cbir_BS %d (%d) %d RW-SW)" % (len(pathway_sets['Gcbir_bs']) - len(pathway_sets['GSWRW-cbir']), len(pathway_sets['GSWRW-cbir']), len(pathway_sets['GSW-RW']) - len(pathway_sets['GSWRW-cbir']))

    print "-" * 50
    print "pathway sets:"
    pathway_sets['allCbir']   = pathway_sets['cbir_meth'] & pathway_sets['cbir_bs']
    pathway_sets['pol-cbir']  = pathway_sets['Polistes']  & pathway_sets['cbir_bs']
    pathway_sets['RWQ-cbir']  = pathway_sets['RW-Q']  & pathway_sets['cbir_bs']
    pathway_sets['SWQ-cbir']  = pathway_sets['SW-Q']  & pathway_sets['cbir_bs']
    pathway_sets['SWRW-cbir'] = pathway_sets['SW-RW'] & pathway_sets['cbir_bs']
    pathway_sets['pol_QW-cbir'] = pathway_sets['Pol_QW'] & pathway_sets['cbir_bs']

    print "(Cbir_BS %d (%d) %d Cbir_meth)" % (len(pathway_sets['cbir_bs']) - len(pathway_sets['allCbir']), len(pathway_sets['allCbir']), len(pathway_sets['cbir_meth']) - len(pathway_sets['allCbir']))
    print "(Cbir_BS %d (%d) %d Polistes)" % (len(pathway_sets['cbir_bs']) - len(pathway_sets['pol-cbir']), len(pathway_sets['pol-cbir']), len(pathway_sets['Polistes']) - len(pathway_sets['pol-cbir']))
    print "(Cbir_BS %d (%d) %d RW-Q)" % (len(pathway_sets['cbir_bs']) - len(pathway_sets['RWQ-cbir']), len(pathway_sets['RWQ-cbir']), len(pathway_sets['RW-Q']) - len(pathway_sets['RWQ-cbir']))
    print "(Cbir_BS %d (%d) %d SW-Q)" % (len(pathway_sets['cbir_bs']) - len(pathway_sets['SWQ-cbir']), len(pathway_sets['SWQ-cbir']), len(pathway_sets['SW-Q']) - len(pathway_sets['SWQ-cbir']))
    print "(Cbir_BS %d (%d) %d RW-SW)" % (len(pathway_sets['cbir_bs']) - len(pathway_sets['SWRW-cbir']), len(pathway_sets['SWRW-cbir']), len(pathway_sets['SW-RW']) - len(pathway_sets['SWRW-cbir']))

    Gpol = pathway_sets['Gpol_QW-cbir']
    Gcbir = pathway_sets['Gcbir_bs']
    GSWQ =  pathway_sets['GSW-Q']
    Pol = pathway_sets['pol_QW-cbir']
    Cbir = pathway_sets['cbir_bs']
    SWQ =  pathway_sets['SW-Q']

    figure, axes = plt.subplots(1, 2)
    venn3([Gpol, Gcbir, GSWQ], ('Polistes', 'Cerapachys', 'Apis'), ax=axes[0])
    venn3([Pol, Cbir, SWQ], ('Polistes', 'Cerapachys', 'Apis'), ax=axes[1])
    plt.show()

########################################################################################

if __name__ == '__main__':

    parser = define_arguments()

    args = parser.parse_args()

    run_arguments(args)