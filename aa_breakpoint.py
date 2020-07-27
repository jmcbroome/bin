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
# 2/6/19 Finalized changes for 2/5, implemented aa comparison algorithm, began testing
# 2/15/19 Took a week off to develop tad2tad, returning to testing. Extremely slow runtime, rework of control generation to use a cs set to cut down on repetitive runs.
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
        self.parser = argparse.ArgumentParser(description = 'Program that takes an axt paired alignment file containing BOS annotations and a bed file containing TADs and returns the significance of proportional enrichment or correlation values.', 
                                             add_help = True,
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -a axt_file -b bed_file [other options*]' 
                                             )
        self.parser.add_argument('-a', '--axt_file', help = "Name of the file containing the Point Of Interest (POI) data in the format Spec/chr/p1/p2/spec2/chr/p1/p2.")
        self.parser.add_argument('-b', '--bed_file', help = "Names of the BED annotation file containing the location of called tad boundaries.")
        self.parser.add_argument('-v', '--verbose', type = bool, default = False, help = "Set to True to print status updates.")
        self.parser.add_argument('-o', '--output_file', help = "Name the file where statistics will be printed. Default is stdout.", default = '')
        self.parser.add_argument('-l', '--chromosome_lengths', help = "File containing the list of lengths of each chromosome. Defaults to Macaque layout.", default = '')
        self.parser.add_argument('-c', '--chain', type = int, help = "Set to a number of monte carlo samples to generate for control distributions. Default 50", default = 50)
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

def axtRead(filename, minq = 0):
    '''
    Read in an axt pairwise alignment file, just retaining the location of BOS's with an alignment score above a quality threshold (default, no threshold)
    '''
    with open(filename) as axt:
        boss = []
        for entry in axt.readlines():
            #axt is formatted similarly to a fasta file. I only want the second through fourth lines of each header and the quality score
            #some entries are headers, some are blank, some are genome sequences
            #skip anything that's not a header.
            if len(entry.strip().split()) > 1: #blank lines are len 0, sequence lines are len 1
                entry = entry.strip().split()
                if entry[2] != entry[5] and entry[3] != entry[6] and int(entry[8]) > minq: #double check that its actually a bos and above quality threshold
                    boss.append([entry[1], int(entry[2]), int(entry[3])])
    return boss 
                
def bedRead(filename):
    '''
    Reads the bed file into a dictionary of dictionaries (first dict key(chrom):value(dict), second dict key(start):value(end))
    '''
    with open(filename) as bed:
        bedDict = {}
        for entry in bed.readlines()[1:]:
            entry = entry.strip().split()
            tmp = bedDict.get(entry[0], {})
            try: #there are some lines that are not actual info, if you run into one, just keep going.
                tmp[int(entry[1])] = int(entry[2])
                bedDict[entry[0]] = tmp
            except:
                continue
    return filename, bedDict

def lenRead(filename):
    with open(filename) as lens:
        return {entry[0]:entry[1] for entry in lens.strip().split()}

class breakPoint:
    def __init__(self, AXTName, AXTarray, BEDdict, out, chromlens, chain = 50, verbose = False):
        '''
        A class that automatically performs a series of operations on a trio of structures representing genomic element stretches sorted by chromosome, individual base points of interest sorted again by chromosome, and lengths of chromosomes on initialization.
        The output is a file containing a statistical summary of the enrichment value of the elements around the points of interest.

        A summary of terms, attributes, and methods:

        Note: Let "coordinate" for this script be the tuple (chrom, base number) where chrom is a string and base number is an integer, representing the chromosome and the location within that chromosome for a base of interest
        
        Input AXTarray is a list of lists containing in each sublist the trio [chr, pos1, pos2] representing a region of BOS vs some other species. 
        Input BEDdict is a dictionary with key of chromosome and value of pairs of boundary bins.
            This does not retain pairing information, which may or may not be important for identification. I think I will leave out pairwise analysis for now but that would be a significant next step for examining, once the different in the inversion types is tested for.
        Input minw and maxw are the minimum and maximum window sizes used to analyze enrichment at each POI. It will always use 10 increments spread evenly between these two points (10 window sizes).
        Input out is a string representing the intended name of the output file
        Input verbose is an optional argument that prints status updates
        
        self.bed is a set containing many tuples of coordinate values of areas of break of synteny.
        self.chromlen is a simple dictionary of chromosome:length (str:int)
        self.controlchain is a dictionary with key chromosome and value of a list of enrichment values pulled from a random distribution of 100 "pois" to act as a statistical control (Monte Carlo)

        self.count is the workhorse of enrichment calculation; for a given window size and poi, it will take all bases within the window range and check for membership in self.gff, and return the proportion that are.
        self.controldist uses self.count to generate a monte carlo set of enrichment values for statistical comparison
        self.test takes a set of values (either distance or proportion) and compares it to the monte carlo distribution for its home chromosome, returning a single pvalue

        Graphs I'll want to make at some point:
        A graph of the distribution of insulator elements across the genome (check if its actually normal or not, assess value of testing function)
        A combined visualization showing inversions bps and insulator elements across each chromosome
        '''
        import math
        if chromlens != '':
            self.chromlen = chromlens
        else:
            self.chromlen = {'chr1':225584828, 'chr2':204787373, 'chr3':185818997, 'chr4':172585720, 'chr5':190429646, 'chr6':180051392, 'chr7':169600520, 'chr8':144306982, 'chr9':129882849, 'chr10':92844088, 'chr11':133663169, 'chr12':125506784, 'chr13':108979918, 'chr14':127894412, 'chr15':111343173, 'chr16':77216781, 'chr17':95684472, 'chr18':70235451, 'chr19':53671032, 'chr20':74971481, 'chrX':149150640, 'chrY':11753682} #values as of rheMac8 assembly
        self.controlchain = {}
        #here I'm writing fresh code to perform the enrichment overlap calcuations. 
        #for each area of interest in the axt
        pvals = {}
        #instead of generating a chainset for each entry in axt, better to create a generic set and round to the nearest when testing.
        cprops = {}
        cs = [50, 500, 1000, 2500, 5000]
        bcs = [10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000] #bcs are huge regions that I shortcut calculate based on cs[-1]
        for csize in cs: #bears more entries if necessary.
            for chro in BEDdict.keys():
                cprops[(chro, csize)] = [self.rprop([0, csize], BEDdict[chro], chro) for _ in range(chain)]
        for bcsize in bcs:
            cprops[(chro, bcsize)] = [math.trunc(prop * (bcsize/5000)) for prop in cprops[(chro, 5000)]] #scaling up 5000 entry for bigger region guesses. 
        for entry in AXTarray:
            try: #adding a try statement to skip over entries in the axt that correspond to chromosomes not covered by the bed
                prop = self.pcount(entry[1:], BEDdict[entry[0]])
                #now I need to generate monte carlo chain values.
                if (entry[2]-entry[1]) > bcs[-1] * 2: #if its too big for reasonable comparison...
                    print("Too large; add entry to control. {}".format((entry[2]-entry[1]))) #if I get many of these errors, I will see an overrepresentation of these guys.
                    #hell I might see all sorts of artifact pvals. if I see way too many significant things I'm going to have to consider how to take into account the difference, or fine-grain the controls more.
                cprop = cprops[(entry[0], min(cs + bcs, key=lambda x:abs(x-entry[2]+entry[1])))]
                pvals[(entry[0], (entry[1], entry[2]))] = self.test(cprop, prop)
            except KeyError:
                pass
        #then to printing.
        with open(out + '_' + AXTName + '.txt', 'w+') as outf:
            print("##################################################", file = outf)
            for chro in sorted(list(set([key[0] for key in list(pvals.keys())]))):
                for startend in sorted(list(set([key[1] for key in list(pvals.keys())]))):
                    try:
                        print("Chromosome: {}\tRegionSize: {}\t Pvalue: {}".format(chro, startend, pvals[(chro, startend)]), file = outf)
                    except:
                        continue #some stuff doesn't want to print right. Just skip it if its a nogo.
        print("Done.")
   
    def pcount(self, axt_region, bed_dict):
        '''
        Function that implements the algorithm descirbed in the brainstorm document for a-a overlap comparison. 
        Axt_region should be a simple pair of points in a list, bed_dict should be a dictionary representing the subset of regions local to the chromosome of the axt pair.
        '''
        check_end = True
        positive = 0
        base = axt_region[0]
        ends = set(bed_dict.values()) #save these in a set for rapid checking if I've hit an end before I've reached any starts.
        while True:
            if check_end == True:
                if base in bed_dict:
                    if bed_dict[base] > axt_region[1]: #if this is the first region and it goes beyond the end
                        positive += axt_region[1] - base
                        break
                    else: 
                        check_end = False
                        positive += bed_dict[base] - base
                        base = bed_dict[base]
                elif base in ends: #if you run into the end of a region that started before the test area. Assumption of no overlap
                    check_end = False
                    positive += base - axt_region[0]
                    base += 1 #count it but keeping moving incremental
                else: #if there's nothing, just keep moving
                    base += 1
                    if base >= axt_region[1]: #unless you've left without hitting anything ofc.
                        break
            else:
                if base in bed_dict:
                    if bed_dict[base] > axt_region[1]:
                        positive += axt_region[1] - base
                        break
                    else:
                        positive += bed_dict[base] - base
                        base = bed_dict[base]
                else:
                    base += 1
                    if base >= axt_region[1]:
                        break
        return positive

    def rprop(self, axt_region, bed_dict, chro):
        '''
        Function that generates a random size-equivalent to the axt region anywhere in the same chromosome and runs pcount to get control proportion values.
        '''
        import random
        rlen = axt_region[1] - axt_region[0]
        rstart = random.randint(10000, self.chromlen[chro] - 10000 - rlen) #assumption that chromosome is longer than 10000 bp, removes the telomeres from consideration ala carbone
        return self.pcount([rstart, rstart + rlen], bed_dict)

    def test(self, cdist, envals):
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
        print("Reading AXT annotation.")
    axt = axtRead(myCommandLine.args.axt_file)
    if myCommandLine.args.verbose:
        print("Reading BED file.")
    bname, bed = bedRead(myCommandLine.args.bed_file)
    if myCommandLine.args.chromosome_lengths != '':
        lens = lenRead(myCommandLine.args.chromosome_lengths)
    else:
        lens = ''
    breakPoint(bname, axt, bed, myCommandLine.args.output_file, lens, myCommandLine.args.chain, myCommandLine.args.verbose)

if __name__ == "__main__":
    main()