#!/usr/bin/env python3

#Jakob McBroome, 11/27/18

#101140 is inversion test file.
#lets define a bin set class

#1/28/19- added mirroring to dense matrix (NxN) writing for compatibility with IC-Finder. Added Windowed option to further subset and do high-resolution contact mapping with reasonable runtimes.
class HiCBin:
    '''
    Binner object that takes filename and minimum and sorts by chromosome on initiation. 
    Returns bin frequencies of a chromosome with self.binner(chromosome, size of bin).
    '''
    def __init__(self, filename, minimum, window = [0,1]):
        #first step is making a dictionary with key chromosome and value contact in that chromosome
        #after this each step occurs per chromosome in the set
        #i guess the alternative would be to pick a chromosome to bin? might be better
        self.chromosomes = {}
        with open(filename) as file:
            raw = file.readlines()
            for entry in raw:
                elist = entry.strip().split()
                if abs(int(elist[2])-int(elist[1])) > minimum and window == [0,1]:
                    if elist[0] not in self.chromosomes.keys():
                        self.chromosomes.update({elist[0]:[(int(elist[1]), int(elist[2]))]})
                    else:
                        self.chromosomes[elist[0]].append((int(elist[1]), int(elist[2])))
                elif abs(int(elist[2])-int(elist[1])) > minimum and window != [0,1]:
                    if int(window[0])<int(elist[1])<int(window[1]) and int(window[0])<int(elist[2])<int(window[1]):
                        if elist[0] not in self.chromosomes.keys():
                            self.chromosomes.update({elist[0]:[(int(elist[1]), int(elist[2]))]})
                        else:
                            self.chromosomes[elist[0]].append((int(elist[1]), int(elist[2])))

    def binnermatrix(self, chrom, size, window = [0,1], mirror = True):
        '''
        The simpler numpy version of the counter.
        '''
        import numpy as np
        import math
        if window == [0,1]:
            binnum = math.trunc(max(self.chromosomes[chrom])[1] / size) + 1
        else:
            binnum = math.trunc((window[1]-window[0])/size) + 1
        # print(math.trunc(max(self.chromosomes[chrom])[1] / size))
        # print(max(self.chromosomes[chrom]))
        # print(binnum)
        counts = np.zeros([binnum, binnum])
        for entry in self.chromosomes[chrom]:
            # print(entry)
            # try:
            counts[math.trunc((entry[0]-window[0])/size)][math.trunc((entry[1]-window[0])/size)] += 1 
            if mirror:
                counts[math.trunc((entry[1]-window[0])/size)][math.trunc((entry[0]-window[0])/size)] += 1 #both the count and its mirror should be incremented for mirroring
            # except IndexError:
                # print("FailInd = {}".format(entry))
        return counts

    def _heatmapper(self, counts, prefix, filename):
        '''
        Called by binnerdict to generate basic Seaborn contact heatmaps. Remember to generate and clear figures each time I use this package to prevent overwriting issues.
        '''
        import matplotlib.pyplot as plt
        #takes a matrix of (coord):val where val is a pval or a contact and generates a heatmap
        import seaborn as sns
        plt.figure()
        ax = sns.heatmap(counts, vmax = 5)
        out = ax.get_figure()
        out.savefig('{}_{}_contacts.png'.format(prefix, filename))
        plt.clf()

def main():    
    import argparse
    import sys
    parser = argparse.ArgumentParser()
    parser.add_argument('--filename', '-f', help = "File name", default = '')
    parser.add_argument('--minimum', '-m', type = int, help = "Minimum contact distance", default = 0)
    parser.add_argument('--chromosome', '-c', type = str, help = "Name of chromosome of interest", default = "X")
    parser.add_argument('--bin_size', '-b', type = int, help = "Size of bins.", default = 40000)
    parser.add_argument('--freq_matrix','-r', type = bool, help = "If False, returns an NxN matrix, otherwise 3xN^2", default = False)
    parser.add_argument('--output', '-o', help = "Name of output file", default = "binner_out.tsv")
    parser.add_argument('--find_point', '-p', help  = "Set to True to go to bin finding mode.", default = False)
    parser.add_argument('--version', '-v', help = "Set to True to print program version.", default = False)
    parser.add_argument('--mirror', '-i', help = "Mirrors contact frequencies across the chromosome line. Defaults to True", default = True)
    parser.add_argument('--window', '-w', nargs = 2, type = int, help = "Pass a pair of coordinate values to look between to generate a subset matrix. Default is whole chromosome", default = [0,1])
    parser.add_argument('--heatmap', '-e', help = "Set to a string to produce a seaborn heatmap of the binned matrix.", default = '')
    args = parser.parse_args()
    assert args.window[1] > args.window[0]
    if args.version:
        print("BinnerScript v1.2")
    obj = HiCBin(args.filename, args.minimum, args.window)
    test = obj.binnermatrix(args.chromosome, args.bin_size, args.window, args.mirror)
    if args.heatmap != '':
        obj._heatmapper(test, prefix = args.heatmap, filename = args.filename)
    if not args.find_point:
        sys.stdout = open(args.output, 'w+')
        if not args.freq_matrix:
            for row in test:
                for value in row:
                    print(int(value), end="\t")
                print('',end='\n')
        else:
            for colind, row in enumerate(test):
                for rowind, value in enumerate(row):
                    print("{}\t{}\t{}\n".format(colind, rowind, value))
    else:
        import math
        try:
            xcoord = input("x =")
            ycoord = input("y =")
            print('{},{}'.format(math.trunc(int(xcoord) / args.bin_size), math.trunc(int(ycoord) / args.bin_size)))
        except:
            print("Error- Coordinates not included or non-integer inputs.")


if __name__ == "__main__":
    main()