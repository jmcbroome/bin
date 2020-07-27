#!/usr/bin/env python3

#script that implements the class developed in tad2tad_dev.ipynb 
#instead of taking simple breakpoint lists, takes 2D.bed files and extracts the breakpoints for comparison.
#developed for the purpose of examining the distribution of differences between four biological replicates, thusly, it will compare all input files pairwise.

class RCMatrix:
    def __init__(self, control, test, verbose = False):
        import numpy as np
        #first step is to take the breakpoint arrays and convert them into tad-identity hash maps
        if verbose:
            print("Building matrices.")
        self.cmat = self.buildMatrix(control, test)
        self.tmat = self.buildMatrix(test, control).transpose()
        #and attributes to convert mat indeces back to spatial coordinates for graphing and potentially severity calculations.
        if verbose:
            print("Arranging spatial coordinates.")
        self.cspat = {ind:tuple(coord) for ind,coord in enumerate(control)}
        self.tspat = {ind:tuple(coord) for ind,coord in enumerate(test)}
#         self.mmat = np.array([
#             [ self.cmat[row][col] * self.tmat[row][col] for col in range(self.cmat.shape[1])
#             ] for row in range(self.cmat.shape[0])
#         ]) #an experimental matrix currently unused.
        if verbose:
            print("Discovering events.")
        self.events = self.classify()
        
    def bp(self, array):
        #takes the 2xN array and converts it into a hash map of points with values of the tad they belong to, not all points but all points in a tad
        bpd = {}
        for entry in array:
            for point in range(entry[0], entry[1]):
                bpd[point] = tuple(entry)
        return bpd
    
    def buildMatrix(self, main, bparray):
        #main is an array ala the inputs to rcmatrix.
        #it will be compared to a bpd which is one of the post self.bp attributes to generate the composition matrix
#         matrix = {}
#         for entry in main:
#             #entry is a single tuple of len = 2.
#             ctd = self._ctDict(entry, bpd)
#             #maybe do a double hash map instead of a traditional matrix? 
#             matrix[tuple(entry)] = ctd
        import numpy as np
        bpd = self.bp(bparray)
        matrix = np.zeros((len(main), len(bparray)))
        for index, entry in enumerate(main):
            ctd = self._ctDict(entry, bpd)
            for subind, ent in enumerate(sorted([tuple(entry) for entry in bparray])):
                matrix[index][subind] = ctd.get(ent, 0)
        return matrix
    
    def _ctDict(self, a, bpd):
        #a is going to be a single point pair while bpd is a self.c/self.t bpd 
        #finds the composition of a of b for any a,b vector pair.
        #how to do this?
        #first get the chunk that you're decomposing
        #then, starting at the start of that chunk, iterate through all points in the chunk and get their values from the matching dict
        #I'm going to load these values into a dictionary with keys equal to anything that has partial identity and return the dict
        #then when I'm using buildMatrix it will initialize with 0 and then fill in values based on this dict
        ctd = {}
        for point in range(a[0], a[1]):
#             print(point)
#             print(bpd.get(point, 'bound'))
            ctd[bpd.get(point, 'bound')] = ctd.get(bpd.get(point, 'bound'), 0) + 1
        ctd = {key:(val/(a[1]-a[0])) for key, val in ctd.items()}
        return ctd #a dictionary with at least one entry ('bound' = boundary element) with values ranging from 0 to 1.
    
    def classify(self):
        #iterate over both matrices and classify events into categories and severity
        #the resulting dictionary can be used to find events of unusual severity and return their locations.
        #it could also be used to investigate association of events.
        #coordinates are row,col e.g. control chunk, test chunk
        events = {'expansion':[], 'reduction':[], 'merge':[], 'division':[], 'insertion':set(), 'deletion':set()} #an additional entry 'shift':[] is planned at some point, but identification is not yet implemented.
        #the six events. All occur in relation to tads (bigger gap = reduction of tad.)
        assert self.cmat.shape == self.tmat.shape
        for row in range(self.cmat.shape[0]-1):
            for col in range(self.tmat.shape[1]-1): #they should have the same shape.
                #to access row, cmat[index]. for column, cmat[:,index]
                #first check for insertions or deletions.
                ccen = self.cmat[row][col]
                tcen = self.tmat[row][col]
                if sum(self.cmat[row]) == 0:
                    events['deletion'].add(tuple(self.cspat[row]))
                    continue
                if sum(self.tmat[:,col]) == 0:
                    events['insertion'].add(tuple(self.tspat[col]))
                    continue
                if 0 < ccen < 1:
                    if sum(self.cmat[row]) < 1:
                        events['reduction'].append((set([self.cspat[row], self.tspat[col]]), 1-sum(self.cmat[row])))
                    if self.cmat[row][col-1] != 0:
                        events['division'].append(((set([self.cspat[row], self.tspat[col]]), set([self.cspat[row], self.tspat[col-1]])), abs((self.cmat[row][col-1] - ccen) * tcen))) #higher = more significant bc the difference is smaller and the breakpoint moves relatively farther
                    if self.cmat[row][col+1] != 0: #if the groups to either side are nonzero
                        events['division'].append(((set([self.cspat[row], self.tspat[col]]), set([self.cspat[row], self.tspat[col+1]])), abs((self.cmat[row][col+1] - ccen) * tcen)))
                if 0 < tcen < 1:
                    if sum(self.tmat[:,col]) < 1:
                        events['expansion'].append((set([self.cspat[row], self.tspat[col]]), 1-sum(self.tmat[:,col]))) #higher = more significant
                    if self.tmat[row-1][col] != 0:
                        events['merge'].append(((set([self.cspat[row], self.tspat[col]]), set([self.cspat[row-1], self.tspat[col]])), abs(self.tmat[row-1][col] - tcen) * ccen))
                    if self.tmat[row+1][col] != 0: #if the groups to either side are nonzero
                        events['merge'].append(((set([self.cspat[row], self.tspat[col]]), set([self.cspat[row+1], self.tspat[col]])), abs(self.tmat[row+1][col] - tcen) * ccen))
        #I have some redundant entries to prune and shift entries to identify.
        #shift events have the same check as removing duplicate merges but across types.
        #a boundary shift has the coordinate of the moved boundary in the test set only
        #the severity value's sign indicates whether its a forward or backward shift and the value is the change relative to the original size
        #setting aside this pruning code block for now because its a rare issue and I have no idea how to integrate it with the coordinate transition.
#         for eind, entry in enumerate(events['merge']):
#             try:
#                 if entry[0][0] == events['merge'][eind+1][0][1] and entry[1] == events['merge'][eind+1][1]: #same event but reciprocal coordinates.
#                     del events['merge'][eind]
# #                 if entry[0][0] == events['division']
#             except:
#                 print("Error removing entry: {}".format(entry))
#         for eind, entry in enumerate(events['division']):
#             try:
#                 if entry[0][0] == events['division'][eind+1][0][1] and entry[1] == events['division'][eind+1][1]: #same event but reciprocal coordinates.
#                     del events['division'][eind]
#             except:
#                 print("Error removing entry: {}".format(entry))
        return events
    
    def graphSevs(self, pref = 'test', line = [0]): #
        #generate severity graphs for each event type with seaborn and save them.
        #pref is prefix options, line is whether to draw a vertical line (like at a predicted inversion point, say...)
        import matplotlib.pyplot as plt
        import seaborn as sns
        import pandas as pd
        for etype, elist in self.events.items():
            if len(elist) > 0 and not isinstance(elist, set):
                if len(elist[0]) > 1:
                    #generate a matrix frame to hold severity data for each type
                    sevs = [e[1] for e in elist]
                    plt.figure()
                    ax = sns.distplot(sevs)
                    out = ax.get_figure()
                    out.savefig('{}_{}_severities.png'.format(pref, etype))
                    plt.clf()
                    #now that the distribution of severities is generated, lets run it along the chromosome line as a scatterplot.
                    #first convert the list of pairs of coordinates to tuples
                    #for ease I'm going to just for every point that's associated with an event + severity graph that.
                    #so any given event will have a point for each bin it affected.
                    #first generate the frame I'll build the plot from.
                    #x = coordinate of event member
                    #y = severity
                    #to set this up, I'm going to decompose my current structure into a dictionary set with single coord: severity for all events in elist
                    #then I'll convert these to vectors and pass to sns.scatterplot
                    vdict = {}
                    for e in elist:
                        for p in e[0]:
                            #unfortunately there are a few different kinds of structures the coordinates are stored in depending on event so this code gets a little messy.
                            if isinstance(p, set):
                                for _ in range(len(p)):
                                    val = p.pop()
                                    if len(val) <= 1:
                                        vdict[val] = e[1]
                                    else:
                                        for i in range(len(val)):
                                            vdict[val[i]] = e[1]
                            elif isinstance(p, tuple):
                                for i in range(len(p)):
                                    vdict[p[i]] = e[1]
                            elif isinstance(p, int):
                                vdict[p] = e[1]
                    #then restructure into a pandas dataframe. 
                    pandver = pd.DataFrame({'coordinate':list(vdict.keys()), 'severity':list(vdict.values())})
                    plt.figure()
                    ax = sns.regplot(x = 'coordinate', y = 'severity', data = pandver)
                    if line != [0]:
                        for val in line:
                            plt.axvline(x=val)
                    out = ax.get_figure()
                    out.savefig('{}_{}_spatial.png'.format(pref, etype))
                    plt.clf()

    def filter(self, etype, enum = 5):
        #get the top parameter enum most severe events of parameter etype and print them.
        top = sorted(self.events[etype], key = lambda x:x[1], reverse = True)[0:enum]
        print("Event Type: " + etype)
        for occurrence in top:
            print("Location: {}\tSeverity: {}".format(list(occurrence[0]), occurrence[1]))
        
def readBed(filename):
    with open(filename, 'r') as bedfile:
        bps = {}
        for entry in bedfile:
            e =  bps.get(entry.strip().split()[0], [])
            e.append([int(e) for e in entry.strip().split()[1:3]])
            bps[entry.strip().split()[0]] = e
    return bps

def main():
    import argparse
    import sys
    parser = argparse.ArgumentParser()
    parser.add_argument('bed_files', nargs = "+", help = "Set of filenames for pairwise testing.")
    parser.add_argument('-v', '--verbose', type = bool, default = False, help = "Set to True to print status updates.")
    parser.add_argument('-o', '--output_file', default = 't2t_test.txt', help = "Set to a name to print filtered differences to. Defaults to t2t_test.txt")
    parser.add_argument('-p', '--prefix', default = 'test', help = 'Set to a prefix for printed graphics. Defaults to test.')
    parser.add_argument('-n', '--filter_depth', default = 5, help = "Set to an integer to print the top X events of each event type for each pairwise file to the output. Defaults to 5.")
    parser.add_argument('-c', '--chromosome', default = '', help = "Set to a name of a chromosome to only check that chromosome. Checks all chromosomes by default.")
    args = parser.parse_args()
    if args.verbose:
        print("Processing bed files.")
    pbeds = [readBed(bed) for bed in args.bed_files] #pbeds = processed beds, each beds is a dictionary with key chromosome and value list of actual tad entries.
    for bind, bed in enumerate(pbeds): 
        obeds = pbeds[:bind] + pbeds[bind+1:] #the set of all files not including the current under iteration
        if args.verbose:
            print("Evaluating {}".format(args.bed_files[bind]))
        for obed in obeds:
            if args.chromosome == '':
                for chro in bed.keys():
                    if args.verbose:
                        print("Checking chromosome {}".format(chro))
                    if chro in obed.keys():
                        obj = RCMatrix(bed[chro], obed[chro], args.verbose)
                    if args.verbose:
                        print("Generating graphs.")
                    obj.graphSevs(args.prefix)
                    with open(args.output_file, 'w+') as out:
                        if args.verbose:
                            print("Filtering.")
                        sys.stdout = out #chaging stdout for convenience with current implementation of filter
                        obj.filter(args.filter_depth)
                        sys.stdout = sys.__stdout__ #reset the stdout.
            else:
                if args.verbose:
                    print("Checking chromosome {} only.".format(args.chromosome))
                obj = RCMatrix(bed[args.chromosome], obed[args.chromosome], args.verbose)
                if args.verbose:
                    print("Generating graphs.")
                obj.graphSevs(args.prefix)
                with open(args.output_file, 'w+') as out:
                    if args.verbose:
                        print("Filtering.")
                    sys.stdout = out #chaging stdout for convenience with current implementation of filter
                    obj.filter(args.filter_depth)
                    sys.stdout = sys.__stdout__ #reset the stdout.

if __name__ == "__main__":
    main()