#!/usr/bin/env python3

#import
import argparse

#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-g', '--gff', help = 'Path to input gff.')
    parser.add_argument('-p', '--plus', help = 'name of output for plus strand CDSs')
    parser.add_argument('-m', '--minus', help = 'name of output for minus strand CDSs')
    
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    #insert code
    with open(args.gff) as inf:
        with open(args.plus,'w+') as plusout:
            with open(args.minus,'w+') as minusout:
                for entry in inf:
                    if entry[0] != "#":
                        spent = entry.strip().split()
                        if spent[2] == "CDS":
                            #get the parent information.
                            pinfo = spent[8].split(';')[1].strip("Parent=gene-")
                            data = '\t'.join(spent[0:8] + ['gene_id', '"' + pinfo + '";',])
                            if spent[6] == '+':
                                print(data, file = plusout)
                            else:
                                print(data, file = minusout)

if __name__ == "__main__":
    main()