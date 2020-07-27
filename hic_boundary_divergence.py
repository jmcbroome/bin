# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 12:59:00 2020

@author: Jakob McBroome
"""
#%%
import pandas as pd
import numpy as np
import seaborn as sns
#%%
with open('hg38_splitbound.fa') as finf:
    seqs = []
    current = None
    for entry in finf:
        if entry[0] == '>':
            #save the current entry if any.
            if current != None:
                seqs.append(current)
            #initialize an entry
            current = [e.strip(":").strip() for e in entry[1:].split(':') if e.strip(":") != ''] + ['']
        else:
            current[-1] += entry.strip()
#%%
with open('hg38_splitmaps.sam')as sinf:
    maps = []
    current = None
    for entry in sinf:
        if entry[0] != '@':
           spent = entry.strip().split()
           names = [e.strip(':') for e in spent[0].split(':') if e.strip(":") != '']
           maps.append(names + spent[1:])
[[se for se in e if len(se) < 25] for e in maps]
#%%
#now relate these two together into a new structure
#I think I might actually use a class framework for this.
class Boundary:
    def __init__(self, samd, fad):
        #attributes will include the primary sequence, the location and name and cigar data, and the alignment sequences
        name, refchro, loc, alnset = self.read_sam(samd)
        self.name = name
        self.chro = refchro
        self.loc = loc
        self.alns = self.match_alns(alnset,fad)
    
    def read_sam(self, samd):
        #samd is in this case a set of entries from maps
        #with shared names but alternative alignments
        alnd = []
        for entry in samd:
            name, chro, loc, flag, altchro, pos, mapq, cigar = entry[:8]
            mapseq = entry[11]
            alnd.append([flag, altchro, pos, mapq, cigar, entry[12:], mapseq])
        return name, chro, loc, alnd
            
    def match_alns(self, alnd, fad):
        #extract entries in the fasta that have names that match this ones.
        my_alns = []
        for entry in fad:
            name, chro, loc, altchro, altmap, altseq = entry
            if name == self.name:
                #find the full data entry in alnd
                for sd in alnd:
                    #print(altchro, altmap, sd[:-1])
                    if altchro == sd[1] and int(altmap.split('-')[0])+1 == int(sd[2]):
                        #extract the appropriate subsequence of mapq.
                        my_alns.append(sd[:-1] + [altseq])
        return my_alns
#%%
#lets see how many unique boundaries got a mapping.
print(len(set([e[0] for e in maps]))) 
#83 boundaries, okay. 
#what's the distribution on how many there are?
sns.distplot([[e[0] for e in maps].count(ev) for ev in set([e[0] for e in maps])])
#%%
#Most only map a handful of times, and these I trust. One maps 25 times, and that one is probably all repeat crap.
solos = [ev for ev in set([e[0] for e in maps]) if [e[0] for e in maps].count(ev) == 1]
solos
#%%
#solos are the ones had that a single mapping subsection.
#get their de values
solode = []
for sn in solos:
    for m in maps:
        if sn == m[0]:
            solode.append(float(m[-2].split(':')[-1]))
solode    
#%%
#get a similar set of values for the multiple mappers, want to see if they're more or less similar or what
multide = []
mults = [ev for ev in set([e[0] for e in maps]) if [e[0] for e in maps].count(ev) > 1]
for sn in mults:
    for m in maps:
        if sn == m[0]:
            try:
                multide.append(float(m[-2].split(':')[-1]))
            except:
                multide.append(float(m[-3].split(':')[-1]))
multide
sns.distplot(multide)
sns.distplot(solode)
#%%
#next step is to connect to the coolers using the cooler API and extract the bins relevant to each boundary
#then calculate DI/Insulation for each boundary and arrange it all
#then this data table will be my preliminary results for my thesis draft.
#extract the map locations, then for each map location extract the X surrounding bins from each cooler.
import cooler