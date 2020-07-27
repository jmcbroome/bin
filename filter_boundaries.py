#!/usr/bin/env python3

#import
import argparse
import sys
import numpy as np
#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-p', '--percentile', type = float, help = 'Set to a value between 0 or 1 to represent the proportion of the entries that should be returned', default = .5)
    parser.add_argument('-m', '--metric', help = 'Define a value you want to filter a boundary gff based on. Options are delta, pvalue, tsep. Default is delta', default = 'delta')
    parser.add_argument('-l', '--less', type = bool, help = 'Set to True to get the set of entries with a LOWER value. Default False', default = False)
    #I could add entries for in and out, but I'd rather just use stdin/stdout for a script of this small scale
    args = parser.parse_args()
    return args

def parse_gff_line(line):
    chro, method, dtype, start, stop, tsep, null, null2, info = line.strip().split() #null and null2 are useless dots
    bidr, deltar, pvaluer, tsepr = info.split(';')
    bid = bidr[3:]
    delta = float(deltar.strip("delta="))
    pvalue = float(pvaluer.strip("pvalue=")) #tsep is already tracked from its own column.
    tsep = float(tsep)
    return delta, pvalue, tsep, line.strip() #return the original line for reprinting later.

def main():
    args = argparser()
    data = []
    for entry in sys.stdin:
        delta, pv, tsep, line = parse_gff_line(entry)
        data.append((delta, pv, tsep, line))
    #now sort the data by the metric of choice.
    if args.metric == 'delta':
        delts = [k[0] for k in data]
        dt = np.percentile(delts, args.percentile)
        if not args.less:
            goodent = [d[3] for d in data if d[0] > dt]
            for nent in goodent:
                print(nent)
        elif args.less:
            goodent = [d[3] for d in data if d[0] < dt]
            for nent in goodent:
                print(nent)
    #elif args.metric == 'tsep':
     #   delts = [k[2] for k in data]

if __name__ == "__main__":
    main()