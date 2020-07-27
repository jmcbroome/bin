#!/usr/bin/env python3

#script to expand all intervals in a provided bed file by a given distance to either side of the breakpoint.

#import
import argparse
import sys
#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-b', '--bed', help = 'path to bed file to expand intervals. default stdin', default = None)
    parser.add_argument('-l', '--length', type = int, help = 'number of bases to expand. default 5k', default = 5000)
    parser.add_argument('-o', '--output', help = 'name of file to output expanded bed to. default stdout', default = None)
    parser.add_argument('-r', '--reduce', type = bool, help = 'set to true to reduce instead of expand. default false', default = False)
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    if args.bed != None:
        inf = open(args.bed)
    else:
        inf = sys.stdin
    if args.output != None:
        outf = open(args.output, 'w+')
    else:
        outf = sys.stdout
    if args.reduce:
        length = -args.length
    else:
        length = args.length
    for entry in inf:
        spent = entry.strip().split()
        nstart = max(int(spent[1])-length,0)
        nend = int(spent[2])+length #may go over end of chromosome. 
        nentry = '\t'.join([spent[0], str(nstart), str(nend)] + spent[3:])
        print(nentry, file = outf)

if __name__ == "__main__":
    main()