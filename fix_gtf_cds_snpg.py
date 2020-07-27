#!/usr/bin/env python3
#SNPGenie demands that all CDS's be perfect and lifted over CDSs are not. This script will fix that by chopping bits off the end.
#There will be "no stop codon" errors because of this, but whatever.

from sys import stdin
geneid = None
lines = []
for entry in stdin:
    data = entry.strip().split()
    if data[9] != geneid:
        if geneid != None:
            #get the combined cds length for all relevant entries.
            trim = sum([int(d[4])-int(d[3])+1 for d in lines]) % 3 #end of the last minus start of the first. First and last may be same entry. Modulo'd. +1 for each entry because indexing or something.
            lines[-1][4] = str(int(lines[-1][4]) - trim) #the amount to trim off.
            #flush data
            for l in lines:
                print('\t'.join(l))
        #start collecting data for a new gene.
        lines = []
        geneid = data[9]
    #otherwise, just add the current data to the growing list.
    lines.append(data)
    geneid = data[9]