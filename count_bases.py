#!/usr/bin/env python3

#import
import argparse
import statistics as st
#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    #add args
    parser.add_argument('-f', '--fastq')
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    #insert code
    record = False
    linecounts = []

    with open(args.fastq) as inf:
        for entry in inf:
            if entry[0] == '@' or entry[0] == '>':
                record = True
            elif record:
                linecounts.append(len(entry.strip()))
                record = False
    print('Total number of bases: {} Total number of reads {} Mean bases per read: {}'.format(sum(linecounts), len(linecounts), st.mean(linecounts)))
if __name__ == "__main__":
    main()