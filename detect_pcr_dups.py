#!/usr/bin/env python3

import sys
import statistics as st

for entry in sys.stdin:
    spent = entry.strip().split()
    alts = [b for b in spent[4] if b in 'ACGTN.']#.replace("N",'')
    quals = spent[5]#.replace("0",'')
    assert len(alts) == len(quals)
    nalts = ''
    nquals = ''
    for i, base in enumerate(alts):
        if int(quals[i]) > 1 and base != 'N':
            nalts += base
            nquals += quals[i]
    #now iterate through the cleaned pileup and look for any pcr errors with sufficient circle depth to be counted as real.
    #for now, lets say the median distance between any two spots along the line has to be above some threshold.
    distances = {}
    for i,base in enumerate(nalts):
        if base != '.':
            if base not in distances:
                last = i
                distances[base] = []
            else:
                distances[base].append(i-last)
                last = i
    filtered = distances.copy()
    for k,v in distances.items():
        if v == []:
            filtered.pop(k)
    if len(filtered) > 1:
        if any([st.median(dv) < len(dv)/10 for dv in filtered.values()]):
            nent = '\t'.join([spent[0],spent[1],spent[2],spent[3], nalts, nquals])
            print(nent)