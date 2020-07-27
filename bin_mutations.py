#!/usr/bin/env python3

#import
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#define functions/classes
class MutationBinner:
    def __init__(self, filepath):
        self.loose_data = self._read_input(filepath)
        self.bin_scores = {}

    def _read_input(self, path):
        data = {}
        with open(path) as inf:
            for entry in inf:
                chro, loc, ref, var = entry.strip().split()
                if chro not in data:
                    data[chro] = []
                data[chro].append((int(loc), ref.upper() + '->' + var.upper())) #treating masked bases as standard for now
        #ensure data is sorted.
        for chro, locs in data.items():
            data[chro] = sorted(locs) #earliest location per chromosome to latest.
        return data

    def _bin_counts(self, size, efilter = None): #efilter allows you to choose a mutation type with a string of arrow notation. ONLY ONE FOR NOW.
        ccounts = {}
        for chro, locs in self.loose_data.items():
            counts = np.histogram([k for k,v in locs if efilter == None or v == efilter])[0]
            ccounts[chro] = counts
        return ccounts

    def _calculate_discrimination(self,count_vector):
        #basically this divides the distribution along its mean and asks how biomodal the partitions are. Higher score is better.
        #we use the bimodal separation index- mean v1 - mean v2 / 2 (var v1 + var v2), wikipedia/Zhang et al 2003.
        mc = np.mean(count_vector)
        high = [c for c in count_vector if c > mc]
        low = [c for c in count_vector if c <= mc] #ones that are exactly equal to mean can be counted as below for this.
        di = (np.mean(high)-np.mean(low))/(2*(np.std(high)-np.std(low)))
        return di

    def find_best_bin(self, start, stop, step, efilter = None, ethresh = 100):
        scores = {}
        for size in np.arange(start, stop, step):
            scores[size] = []
            for chro, cvec in self._bin_counts(size, efilter).items():
                if len(cvec) >= ethresh:
                    scores[size].append(self._calculate_discrimination(cvec))
        #convert size score sets to averages and return the best one.
        self.bin_scores.update(scores) #saving this as an attribute so I can call this function repeatedly and save all the results.

    def save_best_bins(self,outpath):
        best = sorted([(size, np.mean(scores)) for size, scores in self.bin_scores.items()],reverse = True,key = lambda x:x[1])
        size = best[0]
        with open(outpath,'w+') as outf:
            for chro, locs in self.loose_data.items():
                start = 0
                stop = size
                count = 0
                for l,t in locs:
                    if start <= l < stop:
                        count += 1
                    elif l >= stop:
                        print('\t'.join([chro, start, stop, count]),file = outf)
                        count = 0
                        start += size
                        stop += size

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-m', '--mutations', help = 'Path to the mutation "vcf" file (chro, point, ref, alt))')
    parser.add_argument('-o', '--output', help = 'Name for output bedgraph which will contain the best set of bins.')
    parser.add_argument('-b', '--bins', nargs = 3, type = int, help = '3 integer values representing the smallest bin, largest bin, and step size between those. Not inclusive of the largest.', default = [1000, 11000, 1000])
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    binner = MutationBinner(args.mutations)
    binner.find_best_bin(*args.bins)
    binner.save_best_bins(args.output)

if __name__ == "__main__":
    main()