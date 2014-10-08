#!/usr/bin/env python

import hicluster
import common_path
import time
from numpy import mean
import os, re
import itertools
from multiprocessing import  Process


def from_the_top(datafile, sample_header, iteration):
    filename = os.getcwd() + "/" + str(iteration) +  ".txt"
    data_table = os.path.realpath(datafile)
    cluster = hicluster.Cluster(data_table, exportPath=filename, firstrow=True, genes_as_rows=True, \
                gene_metric='correlation', sample_metric='correlation', \
                gene_method='complete', sample_method='complete', \
                color_gradient='seismic')
    cluster.sample_header = sample_header
    meanlist = cluster.normaliseData(center=True, norm_var=True, log_t=True, sample_norm=False)

    t_dict = hicluster.find_degs(cluster, group1='_FL', group2='_SP')
    q_dict = hicluster.p_to_q([v[1] for v in t_dict.values()], display_on=False)
    # report output to file:
    t_list = []
    out_h = open(cluster.exportPath[:-4] + '.t_test.list', 'w')
    out_h.write('%-12s %-7s %-7s\n' % ('Gene','p-value', 'q-value'))
    for gene in t_dict:
        if q_dict[t_dict[gene][1]] <= 0.05:
            out_h.write( "%-12s %7.4f %.4f %.1f\n" % (gene, t_dict[gene][1],  q_dict[t_dict[gene][1]], t_dict[gene][2]))
            t_list.append(gene)
    out_h.close()

    deg_h = open(cluster.exportPath[:-4] + '.t_test.list', 'rb')
    degs = []
    degs_p = []
    degs_m = []
    for line in deg_h:
        if len(line.split()) != 4:
            continue
        else:
            geneid = line.split()[0]
            degs.append(geneid)
            if float(line.split()[3]) >= 0:
                degs_p.append(geneid)
            else:
                degs_m.append(geneid)
    return degs, degs_p, degs_m, filename[:-4] + ".t_test.list"

def compare_results(smpl_bs, smpl_meth, bsfile, methfile, outfile, t0, icnt):
    degs_bs = from_the_top(bsfile, smpl_bs, "B" + str(icnt))
    degs_meth = from_the_top(methfile, smpl_meth, "M" + str(icnt))

    degs_bs_total = set(degs_bs[0])
    degs_bs_plus  = set(degs_bs[1])
    degs_bs_minus = set(degs_bs[2])
    degs_meth_total = set(degs_meth[0])
    degs_meth_plus  = set(degs_meth[1])
    degs_meth_minus = set(degs_meth[2])

    num_bs = len(degs_bs_total)
    num_meth = len(degs_meth_total)
    num_common = len(degs_bs_total & degs_meth_total)

    num_bs_plus = len(degs_bs_plus)
    num_bs_minus = len(degs_bs_minus)
    num_meth_plus = len(degs_meth_plus)
    num_meth_minus = len(degs_meth_minus)

    both_plus = len(degs_bs_plus & degs_meth_plus)
    both_minus = len(degs_bs_minus & degs_meth_minus)
    plus_minus = len(degs_bs_plus & degs_meth_minus)
    minus_plus = len(degs_bs_minus & degs_meth_plus)
    out_h = open(outfile, 'a')
    res_str = "%-5d %-5d %-5d || B+ %-5d B- %-5d M+ %-5d M- %-5d || ++ %-5d -- %-5d +- %-5d -+ %-5d || %s %s\n"
    out_h.write(res_str % (num_bs,  num_meth, num_common,
        num_bs_plus, num_bs_minus, num_meth_plus, num_meth_minus,
        both_plus, both_minus, plus_minus, minus_plus,
        " ".join(list(smpl_bs)), " ".join(list(smpl_meth))
        ))
    out_h.close()
    t1 = time.time()
    taken = t1 - t0
    ave_speed = taken / (icnt+1)
    t_remaining = (7000 - icnt) * ave_speed
    m, s = divmod(t_remaining, 60)
    h, m = divmod(m, 60)
    print("%d iterations complete. Approx %d:%d:%d remaining\r" % (icnt, h, m, s))
    if num_common == 0:
        os.system('rm ' + degs_bs[3])
        os.system('rm ' + degs_meth[3])

def proceed(bs_smpl, meth_smpl, record_dict):
    "Determine whether this pattern of Forg and Stat has been tested yet"
    generic_bs = find_generic(bs_smpl)
    generic_meth = find_generic(meth_smpl)

    if (generic_bs, generic_meth) in record_dict:
        return False
    else:
        # add to dictionary as a completed combo:
        record_dict[(generic_bs, generic_meth)] = True
        return True

def find_generic(sample_list):
    "Take a sample_header and convert it to '_SP' and '_FL' names only"
    generic_header = []
    for sample in sample_list:
        smp_srch = re.search("(_SP)|(_FL)", sample)
        if smp_srch is not None:
            generic_term = smp_srch.group(0)
        else:
            generic_term = "neither"
        generic_header.append(generic_term)

    return tuple(generic_header)

if __name__ == "__main__":
    parser = hicluster.define_arguments()
    t0 = time.time()

    outfile = os.getcwd() + "/results.info"
    print "Writing results to:", outfile
    out_h = open(outfile, 'w')
    out_h.write("BS   Meth   Both   || B+       B-       M+       M-     || B+M+      B-M-      B+M-     B-M+ \n")
    out_h.close()

    bsfile = '/home/antqueen/genomics/experiments/analyses/PRO20140917_broodswap_controls/batch_correction/batch_correction.Sep22_15.00.batch_corrected.tbl'
    methfile = '/home/antqueen/genomics/experiments/analyses/PRO20140923_methylation_rnaseq/batch_correction/batch_correction.Sep25_11.53.batch_corrected.tbl'
    icnt = 1

    filename = ".".join([os.getcwd(), "0", "list"])
    data_table = os.path.realpath(bsfile)
    clusterbs = hicluster.Cluster(data_table, exportPath=filename, firstrow=True, genes_as_rows=True, \
                gene_metric='correlation', sample_metric='correlation', \
                gene_method='complete', sample_method='complete', \
                color_gradient='seismic')

    sample_header_bs = clusterbs.sample_header


    data_table = os.path.realpath(methfile)
    clustermeth = hicluster.Cluster(data_table, exportPath=filename, firstrow=True, genes_as_rows=True, \
                gene_metric='correlation', sample_metric='correlation', \
                gene_method='complete', sample_method='complete', \
                color_gradient='seismic')

    sample_header_meth = clustermeth.sample_header

    print sample_header_bs
    print sample_header_meth

    # compare using actual data:
    name_gen_bs = sample_header_bs
    name_gen_meth = sample_header_meth
    compare_results(name_gen_bs, name_gen_meth, bsfile, methfile, outfile + ".actual.info", t0, icnt)


    #################### Change this to get more speed!!! ###################
    num_cores = 30
    #########################################################################

    # keep a dictionary of searched groups. For each iteration, check if it
    # has already been searched, and move on if it has.
    record_dict = {}


    # mix it up!
    name_gen_bs = itertools.permutations(sample_header_bs)

    for smpl_bs in name_gen_bs:
        name_gen_meth = itertools.permutations(sample_header_meth)
        for smpl_meth in name_gen_meth:



                # analyse this combo (with multithreading):
                procs = {}
                for p in range(num_cores):
                    while proceed(smpl_bs, smpl_meth, record_dict) is False:
                        try:
                            smpl_meth = name_gen_meth.next()
                        except StopIteration:
                            break
                    icnt += 1
                    resultfile = outfile[:-4] + str(p) + ".info"
                    procs[p] = Process(target=compare_results, args=(smpl_bs, smpl_meth, bsfile, methfile, resultfile, t0, icnt))
                    try:
                        smpl_meth = name_gen_meth.next()
                    except StopIteration:
                        break


                for p in procs:
                    procs[p].start()

                for p in procs:
                    procs[p].join()

    # print actual data again at end (to make finding it easier!):
    print sample_header_bs
    print sample_header_meth




