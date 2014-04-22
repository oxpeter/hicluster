#!/usr/bin/env python

import re
import argparse

import hicluster


def count_names(filename):
    "creates a dictionary of words and their counts"
    
    wordich = {}
    file_h = open(filename, 'rb')
    for line in file_h:
        words = re.findall('([A-Za-z]*)',line)
        for word in words:
            if len(word)>=1:
                if word in wordich:
                    wordich[word] += 1
                else:
                    wordich[word] = 1
    file_h.close()
    return wordich
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parses VCF files. To perform allele ratio test, first run on each VCF file seperately.")
    parser.add_argument("text_files", type=str, help="Comma delimited list of text files for clustering")
    parser.add_argument("-o", "--output", type=str, dest="output", default='text_clust.tbl', help="specify out file")
    args = parser.parse_args()

    # create list of filenames from command line:
    filenames = args.text_files.split(',')
    
    # extract count of words in each filename, also create cumulative dictionary of words:
    alldich = {}
    namedich = {}
    for fn in filenames:
        alldich[fn] = count_names(fn)
        namedich.update(alldich)
    
    # create table:
    out_h = open(args.output, 'w')
    out_h.write("words " + " ".join(namedich.keys()) + "\n")
    for txtfile in filenames:
        out_h.write(txtfile + " ")
        for word in namedich:
            if word in alldich[txtfile]:
                out_h.write(alldich[txtfile][word] + " ")
            else:
                out_h.write("0 ")
        out_h.write("\n")
    out_h.close()
    
        
        
    
    