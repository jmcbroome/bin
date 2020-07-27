#!/usr/bin/env python3

#This is the large coordinator script described in the hic_inversion document that does the actual work of finding and predicting enrichment from inversions.
#Jakob M, 12/17/18

class HiCBin:
    '''
    Binner object that takes filename and minimum and sorts by chromosome on initiation. 
    Returns bin frequencies of a chromosome with self.binner(chromosome, size of bin).
    '''
    def __init__(self, filename, minimum):
        #first step is making a dictionary with key chromosome and value contact in that chromosome
        #after this each step occurs per chromosome in the set
        #i guess the alternative would be to pick a chromosome to bin? might be better
        import math
        self.chromosomes = {}
        self.filename = filename
        with open(filename) as file:
            raw = file.readlines()
            for entry in raw:
                elist = entry.strip().split()
                if int(math.trunc(float(elist[2])))-int(math.trunc(float(elist[1]))) > minimum:
                    if elist[0] not in self.chromosomes.keys():
                        self.chromosomes.update({elist[0]:sorted([(int(math.trunc(float(elist[1]))), int(math.trunc(float(elist[2]))))])})
                    else:
                        self.chromosomes[elist[0]].append(sorted((int(math.trunc(float(elist[1]))), int(math.trunc(float(elist[2]))))))

    def binnermatrix(self, chrom, size, verbose = False):
        '''
        Numpy based counter iterates through the pairs file chromosome delineated in the inits and adds counts to bins based on the coordinates of each entry.
        '''
        import numpy as np
        import math
        edgeval = (max([chro[0] for chro in self.chromosomes[chrom]]), max([chro[1] for chro in self.chromosomes[chrom]]))
        binnum = math.trunc(max(edgeval) / size) + 2
        counts = np.zeros([binnum, binnum])
        for entry in self.chromosomes[chrom]:
            try:
                counts[math.trunc(entry[0]/size)][math.trunc(entry[1]/size)] += 1
            except IndexError:
                pass
                #getting a fair amount of these errors when I don't leave extra rows of bins on the end as edges are clipped off or not binned properly.
                # if verbose:
                    # print("FailInd = {}".format(entry))
        return counts
    
    def binnerdict(self, chrom, size, oprint = '', verbose = False):
        '''
        Runs binnermatrix then converts the numpy matrix into the (coord,coord):val format used by the rest of the script. Optionally, generates a basic heatmap with Seaborn showing bin values by coordinates.
        '''
        np_counts = self.binnermatrix(chrom, size, verbose)
        counts = {}
        for rowind, row in enumerate(np_counts):
            for colind, col in enumerate(row):
                counts[(rowind, colind)] = np_counts[rowind][colind]
        if oprint != '':
            self._heatmapper(matdict = counts, prefix = oprint, filename = self.filename)
        return counts

    def rmat(self, ilocs): #takes iloc in terms of absolute chromosome location
        '''
        Takes the predicted inversion site and alters every contact location within the predicted inversion, regardless of significance of the individual points.
        '''
        nfile = {chrom:val for chrom, val in self.chromosomes.items()}
        nchrom = {chrom:[] for chrom in self.chromosomes.keys()}
        for iloc in ilocs:
            irange = iloc[1:]
            assert irange[1] > irange[0]
            for coord in nfile[iloc[0]]:
                if irange[0] < coord[0] < irange[1] and coord[1] > irange[1]: #x is in and y is out
                    newx = irange[0] + (irange[1] - coord[0])
                    nchrom[iloc[0]].append((newx, coord[1]))
                elif coord[0] < irange[0] and irange[0] < coord[1] < irange[1]: #x is out and y is in
                    newy = irange[0] + (irange[1] - coord[1])
                    nchrom[iloc[0]].append((coord[0], newy))
                elif irange[0] < coord[0] < irange[1] and irange[0] < coord[1] < irange[1]:
                    newx = irange[0] + (irange[1] - coord[0])
                    newy = irange[0] + (irange[1] - coord[1])
                    nchrom[iloc[0]].append((newx, newy))
                else:
                    nchrom[iloc[0]].append(coord) #its out altogether, do nothing
            assert len(nfile[iloc[0]]) == len(nchrom[iloc[0]])
            nfile[iloc[0]] = nchrom[iloc[0]]
        return nfile

    def pairmaker(self, nfile, name = "uninv_test.pairs"):
        '''
        Regenerates a .pairs file based on a chromosome matrix and saves it to the local directory.
        '''
        with open(name, 'w') as out_file:
            for chrom, val in nfile.items():
                for s, e in val:
                    print("{}\t{}\t{}".format(chrom, int(s), int(e)), file = out_file)

    def _heatmapper(self, matdict, prefix, filename):
        '''
        Called by binnerdict to generate basic Seaborn contact heatmaps. Remember to generate and clear figures each time I use this package to prevent overwriting issues.
        '''
        import matplotlib.pyplot as plt
        #takes a matrix of (coord):val where val is a pval or a contact and generates a heatmap
        import numpy as np
        import seaborn as sns
        maxx = max([coord[0] for coord in matdict.keys()]) + 1
        maxy = max([coord[1] for coord in matdict.keys()]) + 1
        vals = np.zeros([maxx, maxy])
        for coord, val in matdict.items():
            vals[coord[0]][coord[1]] = val
        plt.figure()
        ax = sns.heatmap(vals, vmax = 10)
        out = ax.get_figure()
        out.savefig('{}_{}_contacts.png'.format(prefix, filename))
        plt.clf()

class Pstats:
    '''
    Class that performs the comparative statistics that yields a pvalue matrix, uses that matrix to generate a set of contiguous regions of significance, then generates distribution statistics about those regions.
    '''
    def __init__(self, filename, matrix, control, out_name, min_dist, sig, print_group = '', verbose = True):
        '''
        Generates a pvalue matrix and assigns each significant pvalue a group number in initiation.
        '''
        #here is where I need to define self.matrix and self.highsigs.
        import math
        self.filename = filename
        self.inversions = []
        #store filename for saving figures
        self.matrix = matrix
        #self.matrix = count matrix.
        if verbose:
            print("Performing statistical tests.")
        pmatrix = self.pmat(control, out_name, min_dist)
        #self.highsigs = coordinate:groupnum matrix.
        self.highsigs = {coord:-1 for coord, value in pmatrix.items() if value > -math.log(sig, 10)}
        if verbose:
            print("{} Points of Significance found".format(len(self.highsigs)))
        #sigset is sorted set of coordinates by pval
        sigset = [key for key, value in sorted(self.highsigs.items(), key = lambda x:x[1], reverse = True)]
        #first group number. Unassigned coords are labeled -1.
        groupnum = 0
        if verbose:
            print("Assigning groups.")
        for coord in sigset:
            if self.highsigs[coord] == -1:
                self._trace(matrix, coord, groupnum, verbose)
                groupnum += 1
        if print_group != '':
            #not currently written as a heatmap function, though I guess it could be? But it would have few discrete values, not best representation. For figures better as labels on a pval heatmap.
            with open(print_group, 'w+') as out:
                maxy = math.trunc(max([coord[1] for coord in matrix.keys()])) + 1
                for x in range(math.trunc(max([coord[0] for coord in matrix.keys()])) + 1):
                    for y in range(maxy):
                        print(self.highsigs.get((x,y), -1), end = '\t', file = out)
                    print('', end = '\n', file = out)
    
    def pmat(self, control, out_name, mindist = 0):
        '''
        Does the work of generating the pvalue matrix. Optionally, generates a pvalue heatmap.
        Applies a fishers exact test on the following 2x2 matrix:
        Count inside the bin on the test    Count inside the bin on the control
        Count outside the bin on the test   Count outside the bin on the control
        Gets a pvalue which is associated with that bin coordinate.
        '''
        import math
        from scipy.stats import fisher_exact
        pvec = {}
        matrixsum = sum(self.matrix.values())
        controlsum = sum(control.values())
        for coord, value in [(coord, value) for coord, value in self.matrix.items() if abs(coord[0]- coord[1]) > mindist]:
            oddsratio, pval = fisher_exact([[value, matrixsum-value],[control.get(coord, 0), controlsum-control.get(coord, 0)]]) 
            try:
                pvec[coord] = -math.log(pval, 10)
            except:
                print("Math Domain: {}".format(pval))
                print("Val:{}".format(value))
                pvec[coord] = .00000000000000000000000000000000000001
        if out_name != '':
            self._heatmapper(pvec, out_name, self.filename)
        return pvec

    def _trace(self, matrix, start, groupnum, prev = [], verbose = False, smear = 25):
        '''
        Recursive function that is the workhorse of the region definer. Functions like so:
        For each bin that is rated significant under some arbitrary given metric:
            Assign it a group number
            Look at the bins within some given distance to that bin (horizontal and vertical only). 
            If they are significant as well:
                Start this process over with that bin as the start, giving it the same group number.
            If none are significant, break, increment the group number, and return to the next unassigned significant bin in the list, beginning again.         
        '''
        #starting at the given coord does a significance trace outward and returns the set
        self.highsigs[start] = groupnum
        #then for each cardinal direction...
        while True:
            newbins = []
            #more is added to newbins each time based on a 'smear' metric
            for ext in range(1, smear):
                newbins.extend([(start[0] + ext, start[1]),(start[0] - ext, start[1]),(start[0], start[1] + ext),(start[0], start[1] - ext)])
            binsadded = False
            for nbin in newbins:
                if self.highsigs.get(nbin, -2) == -1:
                    self._trace(matrix, nbin, groupnum, start, verbose)
                    binsadded = True
            #breaks if no bins are added.
            if binsadded == False:
                break

    def _grouplister(self):
        '''
        Short function that inverts the highsig dictionary, so instead of each coordinate as key with value of their group, its each group as key with value of a list of all coordinates in that group.
        '''
        groupd = {}
        for coord, groupn in self.highsigs.items():
            ngroup = groupd.get(groupn, [])
            ngroup.append(coord)
            groupd[groupn] = ngroup
        return groupd

    def _heteroguesser(self, coords, hard_break = (), verbose = False):
        '''
        Method takes a group (a value from groupd) of highsig coordinates and predicts whether the group came from a homozygous or heterozygous inversion.
        Does this by comparing proportions. In a heterozygotic inversion we have a blobby shape with no particular biases; in a homozygous shape we have the butterfly pattern.
        First guesses at a center using the basic method (difference between smallest x and largest x divided by 2, same for y) then assigns all points in the coord set to a quadrant based on that center
        If the vast majority of points are in quadrants 2 and 4 and not 1 and 3, guess homozygous
        If they're in any other distribution, guess heterozygous
        Takes the naive center prediction and runs
        '''
        if hard_break == ():
            from scipy.optimize import minimize
            xmin = min([coord[0] for coord in coords])
            xmax = max([coord[0] for coord in coords]) 
            ymin = min([coord[1] for coord in coords]) 
            ymax = max([coord[1] for coord in coords])
            ncenter = ((xmax-xmin)/2 + xmin, (ymax-ymin)/2 + ymin)
            xcen = ncenter[0]
            ycen = ncenter[1]
            #this is where I get fancier with the center prediction
            #its going to try to minimize the total distance function based on the group coords and a predicted bp.
            newest = minimize(self.distcomp, ncenter, args = coords, method = "Nelder-Mead", options = {'maxiter':5000,'maxfev':5000})
            # print(newest)
            #     print('Original Guess {}-{}'.format(xcen, ycen))
            #     print('Optimized Guess: {}, {}, {}, {}'.format(newest.fun, newest.nfev, newest.success, newest.x))
            xcen = newest.x[0]
            ycen = newest.x[1]
        else:
            xcen = hard_break[0] 
            ycen = hard_break[1] 
        quadrants = {q:[] for q in [1,2,3,4]}
        for coord in coords:
            if coord[0] < xcen and coord[1] < ycen:
                quadrants[2].append(coord)
            elif coord[0] > xcen and coord[1] < ycen:
                quadrants[1].append(coord)
            elif coord[0] > xcen and coord[1] > ycen:
                quadrants[4].append(coord)
            elif coord[0] < xcen and coord[1] > ycen:
                quadrants[3].append(coord)    
        #replace this line with a test of some kind? 
        if (len(quadrants[2]) + len(quadrants[4])) / (len(quadrants[1]) + len(quadrants[3]) + 1) < 2:
            #then its probably heterozygous
            typec = "Heterozygous"
        else:
            typec = "Homozygous"    
        return (xcen, ycen, typec)

    def distcomp(self, bpoint, coords):
        import math
        dist = 0
        for coord in coords:
            xdist = abs(bpoint[0] - coord[0])
            ydist = abs(bpoint[1] - coord[1])
            #lets try adding a factor of the contact value of the bin so it weights more heavily towards the highest expression bins, see if that helps accuracy
            dist += math.sqrt(xdist ** 2 + ydist ** 2) * self.matrix[coord] ** 6

        return dist

    def _newgroupstat(self, out, chrom, region, binsize = 150000, flipped = False, hard_break = (), verbose = False, minimal = False):
        '''
        Does the printing work, as well as calculating the corners of the region in a quadrilateral fashion. 
        '''
        groupd = self._grouplister()
        with open(out, 'a+') as ofile:
            if not minimal:
                print("Name of File: {}\tChromosome: {}".format(self.filename, chrom), file = ofile)
                if flipped:
                    print("Post-Flipping Analysis", file = ofile)
                print("#######################################", file = ofile)
                for group, coords in groupd.items():
                    pset = [self.matrix[coord] for coord in coords]
                    xmin = min([coord[0] for coord in coords]) * binsize
                    xmax = max([coord[0] for coord in coords]) * binsize
                    ymin = min([coord[1] for coord in coords]) * binsize
                    ymax = max([coord[1] for coord in coords]) * binsize
                    center = self._heteroguesser(coords, hard_break, verbose)
                    if len(coords) > region:
                        print("Group #: {}\tSize: {} Bins\tHighest P: {}\tAvg P: {}".format(group, len(coords), max(pset), sum(pset)/len(coords)), file = ofile)
                        if xmax < ymin * .975:
                            print("S1: {}\tS2: {}\tE1: {}\tE2: {}".format(xmin, xmax, ymin, ymax), file = ofile)
                            print("Predicted Start: {}\tPredicted End: {}\tType: {}".format(center[0] * binsize, center[1] * binsize, center[2]), file = ofile)
                            print('', end = '\n', file = ofile)
                            #adding a line to save the inversion start/ends as an object attribute, accessible for binner inverting.
                            #abusing object attributes but whatever.
                            if center[2] == "Homozygous": #inverting heterzygous inversions is not effective.
                                self.inversions.append([chrom, center[0] * binsize, center[1] * binsize])
            else:
                for group, coords in groupd.items():
                    pset = [self.matrix[coord] for coord in coords]
                    xmin = min([coord[0] for coord in coords]) * binsize
                    xmax = max([coord[0] for coord in coords]) * binsize
                    ymin = min([coord[1] for coord in coords]) * binsize
                    ymax = max([coord[1] for coord in coords]) * binsize
                    center = self._heteroguesser(coords, hard_break, verbose)
                    if len(coords) > region:
                        if xmax < ymin * .975:
                            print(chrom + '\t' + str(center[0]*binsize)[:-2] + '\t' + str(center[1]*binsize)[:-2], file = ofile)
                            if center[2] == "Homozygous":
                                self.inversions.append([chrom, center[0] * binsize, center[1] * binsize])

    def _heatmapper(self, matdict, prefix, filename):
        '''
        Generates the pvalue version of the binner class heatmap. 
        '''
        import matplotlib.pyplot as plt
        import numpy as np
        import seaborn as sns
        import math
        sns.reset_orig()
        maxx = max([coord[0] for coord in matdict.keys()]) + 1
        maxy = max([coord[1] for coord in matdict.keys()]) + 1
        valsa = np.zeros([maxx, maxy])
        for coord, val in matdict.items():
            valsa[coord[0]][coord[1]] = math.trunc(abs(val))
        plt.figure()
        axe = sns.heatmap(valsa, vmax = 10, cmap="Spectral_r")
        oute = axe.get_figure()
        oute.savefig('{}_{}_pvals.png'.format(prefix, filename))
        plt.clf()

def main():
    import argparse, math
    parser = argparse.ArgumentParser()
    parser.add_argument("test_files", nargs = "+", help = "Names of files for testing.")
    parser.add_argument("--control_file", '-f', help = 'Name of file as control.', required = True)
    parser.add_argument("--chromosome", '-r', help = "Set a single chromosome to check only. Default checks all.", default = '')
    parser.add_argument("--bin_size", '-b', type = int, help = "Set the bin size to use for this analysis.", default = 40000)
    parser.add_argument("--minimum_contact", '-c', type = int, help = "Set the minimum distance between contact for consideration.", default = 0)
    parser.add_argument("--minimum_distance", '-d', type = int, help = "Set the minimum distance represented by bins for differential testing", default = 2)
    parser.add_argument("--minimum_significance", '-g', type = float, help = "Set to a minimum pvalue for consideration as significant- 1x10^-6 is default.", default = .0000001)
    parser.add_argument("--minimum_size", '-s', type = int, help = "Set the minimum size of a contiguous region for representation in the output.", default = 3)
    parser.add_argument("--stats_file", '-o', help = "Name of the output file where the groups statistics will be written", required = True)
    parser.add_argument("--print_contact", '-x', help = "Set to a prefix to print the contact matrix for each test to it as a heatmap.", default = '')
    parser.add_argument("--print_control", '-w', help = "Set to a prefix to print the contact matrix for the control to it as a heatmap.", default = '')
    parser.add_argument("--print_pval", '-y', help = "Set to a prefix to print the pval matrix for each test to it as a heatmap.", default = '')
    parser.add_argument("--print_groups", '-z', help = "Set to a name to print the group matrix for each test to it.", default = '')
    parser.add_argument("--verbose", '-v', type = bool, help = "Set to true to print status updates to the console.", default = False)
    parser.add_argument("--fliptest", '-l', type = bool, help = "Set to true to automatically attempt to remove inversions and repeat the analysis for validation.", default = False)
    parser.add_argument("--hard_break", '-i', nargs = '+', type = int, help = "Give a coordinate pair if you know your exact breakpoint to bypass some steps; also used for testing/debugging.", default = ())
    parser.add_argument("--minimal", '-m', type = bool, default = False, help = "Set to True to print minimal information to the stats output file.")
    args = parser.parse_args()
    #now we start the actual procedure.
    control = HiCBin(args.control_file, args.minimum_contact)
    if args.verbose:
        print("Control binned.")
    for file in args.test_files:
        test = HiCBin(file, args.minimum_contact)
        if args.verbose:
            print("{} binned.".format(file))
        if args.chromosome == '':
            chr_test = ["X", "2", "3"]
        else:
            chr_test = [args.chromosome]
        for chrom in chr_test:
            if args.verbose:
                print("Testing Chromosome {}".format(chrom))
            control_contacts = control.binnerdict(chrom, args.bin_size, args.print_control, args.verbose)
            test_contacts = test.binnerdict(chrom, args.bin_size, args.print_contact, args.verbose)
            obj = Pstats(file, test_contacts, control_contacts, args.print_pval, args.minimum_distance, args.minimum_significance, args.print_groups, args.verbose)
            if args.hard_break != ():
                obj._newgroupstat(args.stats_file, chrom, args.minimum_size, args.bin_size, hard_break = [int(math.trunc(hbreak/args.bin_size)) for hbreak in args.hard_break], verbose = args.verbose, minimal = args.minimal)
            else:
                obj._newgroupstat(args.stats_file, chrom, args.minimum_size, args.bin_size, verbose = args.verbose, minimal = args.minimal)
            if args.fliptest:
                if args.verbose:
                    print("Flipping and rerunning.")
                    print("Generating pairs file.")
                nfile = test.rmat(obj.inversions)
                test.pairmaker(nfile, 'false_'+file)
                newbin = HiCBin('false_' + file, args.minimum_contact)
                newmat = newbin.binnerdict(chrom, args.bin_size, '', args.verbose)
                test._heatmapper(newmat, args.print_contact + "_flip", test.filename)
                ppval = args.print_pval + "_flip"
                uninv = Pstats(file, newmat, control_contacts, ppval, args.minimum_distance, args.minimum_significance, args.print_groups, args.verbose)
                if args.verbose:
                    print("Running flipped stats.")
                if args.hard_break != ():
                    uninv._newgroupstat('flip_'+args.stats_file, chrom, args.minimum_size, args.bin_size, hard_break = [int(math.trunc(hbreak/args.bin_size)) for hbreak in args.hard_break], flipped = True, verbose = args.verbose, minimal = args.minimal)
                else:
                    uninv._newgroupstat('flip_'+args.stats_file, chrom, args.minimum_size, args.bin_size, flipped = True, verbose = args.verbose, minimal = args.minimal)
            if args.verbose:
                print("Complete.")

if __name__ == "__main__":
    main()