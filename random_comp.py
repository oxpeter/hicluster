#!/usr/bin/env python

import hicluster
import common_path
import time
from numpy import mean
import os

def collect_degs(data_file):

    args = parser.parse_args(['-D', data_file,
                '-T', '--no_clustering', '--display_off',
                '--randomise_samples',
                '-fnu',
                '-o', '/home/antqueen/genomics/experiments/analyses/PRO20140925_brain_expt_comparisons/randomised',
                '-S', '0.05'
                ])


    ## arguments to standard pipeline:
    cluster = hicluster.run_arguments(args)
    #degs = hicluster.make_a_list(cluster.exportPath[:-4] + '.t_test.list')
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
    return degs, degs_p, degs_m

if __name__ == "__main__":
    parser = hicluster.define_arguments()
    t0 = time.time()

    num_bs = []
    num_meth = []
    num_common = []

    num_bs_plus = []
    num_bs_minus = []
    num_meth_plus = []
    num_meth_minus = []

    both_plus  = []
    both_minus = []
    plus_minus = []
    minus_plus = []
    print "Writing results to:", os.getcwd() + "/results.info"
    out_h = open(os.getcwd() + "/results.info", 'w')
    out_h.write("iter BS Meth  Both   || B+       B-       M+       M-    || B+M+      B-M-       B+M-     B-M+ \n")
    out_h.close()
    for i in range(1000):
        degs_bs = collect_degs('/home/antqueen/genomics/experiments/analyses/PRO20140917_broodswap_controls/batch_correction/batch_correction.Sep22_15.00.batch_corrected.tbl')
        degs_meth = collect_degs('/home/antqueen/genomics/experiments/analyses/PRO20140923_methylation_rnaseq/batch_correction/batch_correction.Sep25_11.53.batch_corrected.tbl')

        degs_bs_total = set(degs_bs[0])
        degs_bs_plus  = set(degs_bs[1])
        degs_bs_minus = set(degs_bs[2])
        degs_meth_total = set(degs_meth[0])
        degs_meth_plus  = set(degs_meth[1])
        degs_meth_minus = set(degs_meth[2])

        num_bs.append(len(degs_bs_total))
        num_meth.append(len(degs_meth_total))
        num_common.append(len(degs_bs_total & degs_meth_total))

        num_bs_plus.append(len(degs_bs_plus))
        num_bs_minus.append(len(degs_bs_minus))
        num_meth_plus.append(len(degs_meth_plus))
        num_meth_minus.append(len(degs_meth_minus))

        both_plus.append(len(degs_bs_plus & degs_meth_plus))
        both_minus.append(len(degs_bs_minus & degs_meth_minus))
        plus_minus.append(len(degs_bs_plus & degs_meth_minus))
        minus_plus.append(len(degs_bs_minus & degs_meth_plus))
        out_h = open(os.getcwd() + "/results.info", 'a')
        res_str = "%3d) %-5d %-5d %-5d || B+ %-5d B- %-5d M+ %-5d M- %-5d || ++ %-5d -- %-5d +- %-5d -+ %-5d\n"
        out_h.write(res_str % (i+1, num_bs[-1],  num_meth[-1], num_common[-1],
            num_bs_plus[-1], num_bs_minus[-1], num_meth_plus[-1], num_meth_minus[-1],
            both_plus[-1], both_minus[-1], plus_minus[-1], minus_plus[-1]
            ))
        out_h.close()
        t1 = time.time()
        taken = t1 - t0
        ave_speed = taken / (i+1)
        t_remaining = (1000 - i) * ave_speed
        m, s = divmod(t_remaining, 60)
        h, m = divmod(m, 60)
        print("%d iterations complete. Approx %d:%d:%d remaining\r" % (i+1, h, m, s))
    out_h = open(os.getcwd() + "/results.info", 'a')
    out_h.write("MEAN:\n")
    out_h.write(res_str % (i+1, mean(num_bs),  mean(num_meth), mean(num_common),
            mean(num_bs_plus), mean(num_bs_minus), mean(num_meth_plus), mean(num_meth_minus),
            mean(both_plus), mean(both_minus), mean(plus_minus), mean(minus_plus)
            ))
    out_h.close()


