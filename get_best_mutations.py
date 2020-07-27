#!/usr/bin/env python3

import sys
import statistics as st
import numpy as np

def get_dindex(altstring):
    distances = {}
    for i,base in enumerate(altstring):
        if base != '.':
            if base not in distances:
                last = i
                distances[base] = []
            else:
                distances[base].append(i-last)
                last = i
    dindex = {}
    for k,v in distances.items():
        if len(v)> 0:
            dindex[k] = st.median(v)
    return dindex

def make_random(length = 100, bases_to_use = 'A', num = 5):
    string = list('.' * length)
    for b in bases_to_use:
        locs = np.random.choice(length,num,replace = False)
        for l in locs:
            string[l] = b
    return ''.join(string)

def perm_index(leng = 100, num = 3, pnum = 1000):
    indeces = []
    for p in range(pnum):
        tstr = make_random(length = leng, num = num, bases_to_use='A')
        index = get_dindex(tstr)
        indeces.append(index['A'])
    return np.percentile(indeces,5)

def make_pileup_line(spent, mind = 10, pcr_duplicate_track = {}):
    #convert a stripped and split mpileup line into a fake vcf line, filling in default values.
    #pileup: chrom loc ref depth vector_of_alts vector_of_quals
    #vcf: chrom loc ID(.) ref alt qual filter info(dp=...)
    #minimalist entry looking at the docs is just "DP" in info, qual is ., id is ., filter is PASS. So try those.
    chrom, loc, ref, ndepth, vector_of_alts, vector_of_quals = spent
    #real depth isn't depth, its the length of the vector of alts without other symbols or Ns, because most Ns are introduced by the consensus builder and not in the original reads.
    # pcr_duplicate_track = {} #using dynamic programming to save compute cycles for this qc measure
    depth = len([v for v in vector_of_alts if v in 'ACGT.'])
    if depth >= mind:
        quality_alts = []
        quality_quals = []
        cleaner = [b for b in vector_of_alts if b in 'ACGTN.']
        assert len(cleaner) == len(vector_of_quals)
        for i,b in enumerate(cleaner):
            if b in 'ACGT.':
                #check if the qual value is high enough.
                if int(vector_of_quals[i]) > 1:
                    quality_alts.append(b)
                    quality_quals.append(vector_of_quals[i])
        depth = len(quality_alts) #update depth
        # for b in 'ACGT': #for all possible bases
            # if quality_alts.count(b) > depth/4: #if that base is more than 25% of seen bases at this point
                # quality_alts = [base for base in quality_alts if base != b] #remove it, it's almost certainly a germline mutation.
        #if I'm removing germline I can apply an additional filter which should remove PCR duplicates
        skip = '' #record no more than one of the bases that will be included here because of pcr duplicate inflation.
        dindeces = get_dindex(quality_alts)
        pcrc = [] #bases which appear to be in a pcr duplicate cluster get a copy here and then skipped by the main iterator
        for base in 'ACGT':
            basecount = quality_alts.count(base)
            if 2 <= basecount <= depth/4: #doesn't make sense to calculate for singletons, which I intend to skip by default now.
                key = (len(quality_alts), basecount)
                if key not in pcr_duplicate_track:
                    pcr_duplicate_track[key] = perm_index(leng = key[0], num = key[1])
                thresh = pcr_duplicate_track[key]
                if dindeces[base] < thresh: #less than 5% chance of getting a cluster like this. Lock this one to 1 instance
                    #print("QC: Base is skipped for clustering")
                    skip += base
                    #still count it once though
                    pcrc.append(base)
                #if a base exists at appropriate levels but stays above the cluster threshold, it's collected in the fixed alts statement below
            else:
                skip += base
        #collect individual bases not counted as a pcr cluster, plus a singleton of each cluster base
        #rebuild a cleaned pileup line 
        capped = []
        final_alts = []
        final_quals = []
        for i, b in enumerate(quality_alts):
            q = quality_quals[i]
            if b not in skip:
                if b in pcrc and b not in capped:
                    capped.append(b) #count one of them but not again
                    final_alts.append(b)
                    final_quals.append(str(q))
                elif b not in pcrc: #not in a cluster, record.
                    final_alts.append(b)
                    final_quals.append(str(q))
            else:
                final_alts.append('.') #ignore these
                final_quals.append(q)
        #the assumption is also that germline variants have been cleaned out by pilon.
        if depth == 0 or len(final_alts) == 0 or len(final_alts) == final_alts.count("."): #nothing but Ns or reference here.
            return None, pcr_duplicate_track
        else:
            pileup = '\t'.join([chrom, str(loc), ref, str(depth), ''.join(final_alts), ''.join(final_quals)])
            return pileup, pcr_duplicate_track
    else:
        return None, pcr_duplicate_track

pdt = {}
for entry in sys.stdin:
    spent = entry.strip().split()
    line, pdt = make_pileup_line(spent, mind = 100, pcr_duplicate_track = pdt)
    if line != None:
        print(line)