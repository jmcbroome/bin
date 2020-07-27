#!/usr/bin/env python3

#import
import argparse
import sys
#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    #add args
    parser.add_argument('-e', '--errors', help = 'path to input file. default is stdin', default = sys.stdin)
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    #insert code
    sumerrors = {}
    for a in "ACGT":
        for b in "ACGT":
            if a != b:
                sumerrors[(a,b)] = 0
    with open(args.errors) as ef:
        for entry in ef:
            spent = entry.strip().split()
            if len(spent) > 4:
                ref = spent[2]
                cir = spent[4]
                if cir in 'ACGT' and ref in 'ACGT':
                    if cir != ref:
                        sumerrors[(ref, cir)] += 1
    print("Error Counts")
    for k,v in sumerrors.items():
        print(k, '\t', v)

if __name__ == "__main__":
    main()