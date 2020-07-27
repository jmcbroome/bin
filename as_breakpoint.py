#!/usr/bin/env python3

########################################################################
# File: aa_breakpoint.py
# Executable: aa_breakpoint.py
# Purpose: Takes a set of genomic features and breakpoints of interest and analyzes distribution and enrichment in different window sizes around those points. 
# Alternately, performs a similar test looking for distances to the nearest instance of the elements for each point of interest.
# stderr: errors and status
# stdout: write to file
#
# Author: Jakob McBroome
########################################################################

########################################################################
# History/Changelog: 
# 1/25/19 Outlined command line structure, wrote read-in functions, wrote __init__() and main()
# 1/26/19 Wrote methods to count proportions and generate random control distributions, added argument to feed in chromosome lengths
# 1/28/19 Restructured how enrichment values are grouped, implemented 2 tailed t test with unequal variance and basic print function
# 1/31/19 Added a 'nearest' mode that seeks the nearest point noted in the elements file and returns the minimum distance (either up or downstream) for analysis of POI against TAD breakpoints or similar purposes
#           This mode includes a simple print statement, its own reader for a two-column TAD boundary set, and an alternate monte carlo control maker.
# 2/1/19 Restructured print statement for proportional enrichment, self.outwrite is deprecated (printing handled in __init__)
# 2/4/19 Altered gff input to take any number of gff files simultaneously, edit print to produce separate output files for each input gff
# 2/5/19 Forked nearest distance code into ss_breakpoint, retained proportional enrichment for aa comparison.
########################################################################

########################################################################
# CommandLine
########################################################################

class CommandLine() :

    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        
        Implements a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program that takes a fasta file and returns a file containing a list of the least represented k-mers and their associated z-scores.', 
                                             add_help = True,
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -f input_file_path [other options*]' 
                                             )
        self.parser.add_argument('-p', '--poi_file', help = "Name of the file containing the Point Of Interest (POI) data in the format ChrX/tF1. If there is a paired point, two more columns /tChrY/tF2, if not, fill with '.'. ")
        self.parser.add_argument('-f', '--gff_file', nargs = '+', help = "Names of GFF3 annotation files containing the location of elements you're interested in the enrichment of.")
        self.parser.add_argument('-v', '--verbose', type = bool, default = False, help = "Set to True to print status updates.")
        self.parser.add_argument('-m', '--minimum_window', type = int, default = 500, help = 'Set a minimum window size value for enrichment evaluation')
        self.parser.add_argument('-n', '--maximum_window', type = int, default = 1200, help = "Set a maximum window size value for enrichment evaluation")
        self.parser.add_argument('-o', '--output_file', help = "Name the file where statistics will be printed. Default is stdout.", default = '')
        self.parser.add_argument('-l', '--chromosome_lengths', help = "File containing the list of lengths of each chromosome. Defaults to Drosophila Melanogaster layout.", default = '')
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

def GFFRead(filename):
    '''
    Function that translates GFF3 formatted input files into an array format, with a row for each base.
    '''
    array = []
    with open(filename) as gff:
        for entry in gff.readlines()[2:]:
            ent = entry.strip().split()
            fixentry = [ent[0], ent[3], ent[4], ' '.join(ent[8:10])]
            array.append(fixentry)
    return [filename[:-4], array] #return the original name of the file sans the .gff extension

def POIRead(filename):
    '''
    Function that translates a file of points of interest into a dictionary format.
    '''
    poidict = {}
    with open(filename) as poi:
        for entry in poi.readlines():
            nent = entry.strip().split()
            poidict[nent[0]] = poidict.get(nent[0], [])  + nent[1:]
    return poidict 

def lenRead(filename):
    with open(filename) as lens:
        return {entry[0]:entry[1] for entry in lens.strip().split()}

def TADRead(filename, binsize = 1):
    '''
    Reads in a simplified TAD breakpoint matrix (2 column, space/tab delineated, all integers). Argument binsize converts bin coordinates to base coordinates as it reads.
    Returns an array of base values representing outer edges of chromatin domains.
    '''
    tads = []
    with open(filename, 'r') as tadbp:
        for entry in tadbp.readlines():
            tads.append(entry[0] * binsize)
            tads.append(entry[1] * binsize)
    return tads

class breakPoint:
    def __init__(self, GFFName, GFFarray, POIdict, minw, maxw, out, chromlens, verbose = False):
        '''
        A class that automatically performs a series of operations on a trio of structures representing genomic element stretches sorted by chromosome, individual base points of interest sorted again by chromosome, and lengths of chromosomes on initialization.
        The output is a file containing a statistical summary of the enrichment value of the elements around the points of interest.

        A summary of terms, attributes, and methods:

        Note: Let "coordinate" for this script be the tuple (chrom, base number) where chrom is a string and base number is an integer, representing the chromosome and the location within that chromosome for a base of interest
        
        Input GFFArray is a list of lists containing the quartet of [chrom, base number1, base number 2 (10 more than 1 for insulators), element type] for every row (element) in the original file
        Input POIdict is a dictionary with key of chromosome and value of all base numbers of interest on that chromosome in a list
            This does not retain pairing information, which may or may not be important for identification. I think I will leave out pairwise analysis for now but that would be a significant next step for examining, once the different in the inversion types is tested for.
        Input minw and maxw are the minimum and maximum window sizes used to analyze enrichment at each POI. It will always use 10 increments spread evenly between these two points (10 window sizes).
        Input out is a string representing the intended name of the output file
        Input verbose is an optional argument that prints status updates
        
        self.gff is a dictionary that uses key coordinate and value insulator type. It contains a key and value for every base that's within an insulator element (at standard size 10, 10 items with consecutive keys represent a single biological element)
        self.chromlen is a simple dictionary of chromosome:length (str:int)
        self.controlchain is a dictionary with key chromosome and value of a list of enrichment values pulled from a random distribution of 100 "pois" to act as a statistical control (Monte Carlo)

        self.gffProcess divides the GFF dictionary of bases by chromosome value for easier accession as well as splitting ranges from the original array into sets of individual bases
            the purpose of this is so that to check whether any base anywhere is an insulator is simply to see if it is a member of the set of keys in self.gff (an O(1) operation)
        self.count is the workhorse of enrichment calculation; for a given window size and poi, it will take all bases within the window range and check for membership in self.gff, and return the proportion that are.
        self.controldist uses self.count to generate a monte carlo set of enrichment values for statistical comparison
        self.test takes a set of values (either distance or proportion) and compares it to the monte carlo distribution for its home chromosome, returning a single pvalue

        Graphs I'll want to make at some point:
        A graph of the distribution of insulator elements across the genome (check if its actually normal or not, assess value of testing function)
        A combined visualization showing inversions bps and insulator elements across each chromosome
        '''
        import math, statistics
        if verbose:
            print("Processing GFF")
        self.gffProcess(GFFarray)
        if chromlens != '':
            self.chromlen = chromlens
        else:
            self.chromlen = {'chr2L':23515712, 'chr2R':25288936, 'chr3L':25288936, 'chr3R':32081331, 'chrX':23544271} #values as of DM6 assembly, chr 2/3/X only for first testing
        self.controlchain = {}
        #for every window size in say 10 increments (might be a better way to do this) starting with the smallest
        if verbose:
            print("Iterating through window sizes.")
        pvals = {}
        for win in range(minw, maxw, math.trunc((maxw-minw)/10)):
            # print(win)
            #for every point of interest in POI        
            if verbose:
                print("Calculating enrichment values with window size {}".format(win))
            for chro, poil in POIdict.items():
                cdist = self.controldist(chro, win)  #a set of values 0-1 generated by equivalent monte carlo chains
                envals = []
                for poi in poil:
                    #count the number of GFFarray points within that region- test every base in the window against membership in gff- get back as proportion.
                    envals.append(self.count(chro, poi, win)) #a value from 0-1, probably near 0
                    ###commenting out some code for handling paired data, since my first iteration is going to treat them as individual points.
                    # if len(poi) == 2: #since my initial set here is inversions I'm expecting paired data but later it might be single points.
                    #     enval1 = self.count(chro, poi[0], win)
                    #     enval2 = self.count(chro, poi[1], win)
                    #     self.POIenrich[(chro, poi)] = [enval1, enval2]
                ###adding this replcaement line as a test for window effects, remove for legit use
                # envals = self.controldist(chro, win, chain = len(envals))

                if verbose:
                    print("Calculating significance for enrichment values with window size {}".format(win))
                #now apply statistical tests to each enrichment value against the control set
                pvals[(chro, win)] = self.test(cdist, envals)
                if pvals[(chro, win)][1] < .05:
                    with open(out + '.txt', 'a+') as outf:
                        print('Chr:{}\tWin:{}\tMeanMC:{}\tEnrich:{}'.format(chro, win, statistics.mean(cdist), envals), file = outf)
        # print(pvals)
        #now prints the name of the gff as part of the output text.
        with open(out + '_' + GFFName + '.txt', 'a+') as outf:
            print("##################################################", file = outf)
            for chro in sorted(list(set([key[0] for key in list(pvals.keys())]))):
                for win in sorted(list(set([key[1] for key in list(pvals.keys())]))):
                        print("Chromosome: {}\tWindow: {}\t Pvalue: {}".format(chro, win, pvals[(chro, win)]), file = outf)
        print("Done.")
   
    def gffProcess(self, gff, eclass = ''):
        '''
        Breaks apart the gff array into subsets by chromosome to allow quicker subsetting and indexing. 
        '''
        #time to get smart with my data structure, if I want this to take less than 10 million years to run.
        #have to set a self.maxgff attribute for control chain generation btw
        #current plan- create a dictionary (hash map) with keys of (chromosome, base) for all bases in any insulator and value of their type
        #then to check if part of a window is an insulator element just see if it returns anything from the hash dict or not
        if len(gff[0]) == 1: #if its a simple column of elements for nearest mode, just retain that list
            self.gff = gff
        else: #if its fancier, make the dict
            self.gff = {}
            for element in gff:
                for base in range(int(element[1]), int(element[2]) + 1):
                    if eclass == '' or element[3] == eclass:
                        self.gff[(element[0], base)] = element[3]
                    
    def count(self, chro, poi, win):
        '''
        Subsets self.gff to a given window size around a given point and calculates the proportion of the total window that is annotated in the GFF.
        '''
        import math
        elems = 0
        nonelems = 0
        for checkbase in range(math.trunc(int(poi)-win/2), math.trunc(int(poi)+win/2)+1):
            if self.gff.get((chro, checkbase), 0) != 0:
                elems += 1
            else:
                nonelems += 1
        return elems/nonelems

    def controldist(self, chro, win, chain = 1000):
        '''
        Calculates a random set of proportions for statistical comparison.
        '''
        import math, random
        #set of enrichment proportions for random points on the chromosome
        #this is the distribution that you get the pvalue by comparing to
        return [self.count(chro, rpoi, win) for rpoi in [math.trunc(random.random() * self.chromlen[chro]) for _ in range(chain)]]

    def test(self, cdist, envals, nearest = False):
        '''
        Performs statistical testing on a set of envals against the control distribution for a given win (assuming same chromosome) and returns a pvalue.
        Should work equally well for proportion testing (all values 0-1) and nearest testing (values are 0 < integers < chromosome length)
        '''
        from scipy.stats import ttest_ind
        return ttest_ind(cdist, envals, equal_var=False)
  
def main(myCommandLine=None):
    '''
    Main() function for inversion analysis, reads in a GFF file, a POI file, and a chromosome length file, then applies breakpoint analysis and prints the results to a file.
    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(myCommandLine)
    #here's where it happens.
    if myCommandLine.args.verbose:
        print("Reading GFF annotations.")
    gffs = [GFFRead(file) for file in myCommandLine.args.gff_file]
    if myCommandLine.args.verbose:
        print("Reading POI file.")
    poi = POIRead(myCommandLine.args.poi_file)
    if myCommandLine.args.chromosome_lengths != '':
        lens = lenRead(myCommandLine.args.chromosome_lengths)
    else:
        lens = ''
    for gname, gfile in gffs:
        breakPoint(gname, gfile, poi, myCommandLine.args.minimum_window, myCommandLine.args.maximum_window, myCommandLine.args.output_file, lens, myCommandLine.args.verbose)

if __name__ == "__main__":
    main()