#!/usr/bin/env python3
#This script takes a pileup and a set of inferred betabinomial parameters
#it infers the body frequency of a mutation at that site (may be an invisible/undetected 'mutation')
#then it calculates the binomial pmf of getting exactly one detection for the given depth of a site
#these values are tabulated by reference, alternative, and context
#and can be summarized to estimate the expected number of singleton detections given the model parameters.

#import
import argparse
import itertools
from scipy.stats import binom
from Bio import SeqIO as sqio

#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-p', '--pileup', help = 'Path to input pileup file.')
    parser.add_argument('-f', '--reference', help = 'Path to reference fasta')
    parser.add_argument('-a', '--alpha', type = float, help = "alpha parameters for the betabinomial.")    
    parser.add_argument('-b', '--beta', type = float, help = 'beta parameter for the betabinomial.')
    parser.add_argument('-o', '--output', help = 'Name of the output table holding counts.')
    parser.add_argument('-c', '--count', type = int, help = 'Set to the number of times seen you want to model for. Default 1, output is annotated with this assumption, change for debugging purposes.', default = 1)
    args = parser.parse_args()
    return args

def calc_freq(a,b,depth,seen):
    return (seen + a) / (depth + a + b)

def get_context(genome, chro, loc):
    try:
        prior = genome[chro][int(loc)-2] #has to do with differences in indexing.
        latter = genome[chro][int(loc)]
    except IndexError:
        #print("Can't access location", chro, loc)
        prior = None
        latter = None
    return prior, latter

def collect_values(args):
    rs = 0 #a tracker variable for the actual number of singletons as a point of comparison.
    pileup = args.pileup
    vd = {} #build the key set of priors, latters, references, alternatives that are possible. Default values to 0.
    for trio in set(itertools.permutations("AAACCCGGGTTT", 3)):
        vd.update({tuple(list(trio) + [a]):0 for a in 'ACGT' if a != trio[1]})
    #get the genome file for context collection.
    genome = sqio.to_dict(sqio.parse(args.reference, format = 'fasta'))
    #iterate through the pileup.
    with open(pileup) as inf:
        for entry in inf:
            chro, loc, ref, depth, alts, quals = entry.strip().split()
            depth = int(depth)
            if depth == 0:
                continue
            alts = [a for a in alts if a in "ACGTN."]
            assert len(alts) == depth
            #get the prior and latter bases.
            prior, latter = get_context(genome, chro, loc)
            if any([b == 'N' for b in [prior, latter, ref]]) or prior == None or latter == None:
                continue
            for a in set(alts):
                if a == 'N' or a == '.':
                    continue
                if alts.count(a) == args.count:
                    rs += 1
                #calculate the binomial pmf of exactly once, using depth and inferred frequency.
                freq = calc_freq(args.alpha,args.beta,depth,alts.count(a))
                #get the binomial pmf for exactly one occurrence.
                prob = binom.pmf(args.count, depth, freq) #may be a very small number.
                #add it to the vd counter.
                vd[(prior, ref, latter, a)] += prob
    return rs, vd

def print_table(vd, outf):
    with open(outf,'w+') as of:
        print('\t'.join(['Prior', 'Reference', "Latter", "Alternative", 'TotalSingleton']), file = of)
        for k,v in vd.items():
            print('\t'.join(list(k)) + '\t' + str(v), file = of)
        
def main():
    args = argparser()
    rs, vd = collect_values(args)
    print_table(vd, args.output)
    #summarize the values and compare to the real singleton tracker.
    #the number of errors should be equal to real number - expectations
    exp = sum(vd.values())
    print("Expected number of singletons:", exp)
    print("Actual number of singletons:", rs)
    prop = (rs - exp) / rs
    print("Proportion of singleton that is error:", prop)

if __name__ == "__main__":
    main()