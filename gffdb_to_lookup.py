#!/usr/bin/env python3
#This hardcoded script looks for a gff database and calculates a custom lookup file
#this file contains, for each chro/loc that is coding and could have a valid frame/cds identified
#the set of mutations which are synonymous and nonsynonymous for that site
#this custom format will be parsed in by scripts down the line, which will use it to divide a pileup into synonymous and nonsynonymous sites.

import gffutils
db = gffutils.FeatureDB('Dsim.db')
translate = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y','TAA':'Stop','TAG':'Stop','TGT':'C','TGC':'C','TGA':'Stop','TGG':'W',
                'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
print('\t'.join(['Chro', 'Loc', 'Ref', 'SynAlts', 'NonAlts']))
def get_synonymity_info(scv):
    syndata = {}
    for i,cb in enumerate(scv):
        syndata[cb] = {k:[] for k in ['syn','non']}
        #first calculate whether this position is first, second, or third
        frame = i%3
        #get the codon this matches appropriately.
        codon = ''.join([scb[1] for scb in scv[i-frame:i+3-frame]]) 
        try:
            refaa = translate[codon]
        except:
            #many codons will contain Ns, and should be skipped.
            #print("Cant translate codon", codon)
            #print(i, frame, i-frame, i+2-frame, i+2-frame - i-frame)
            continue
        #perturb the codon in all ways relevant to this position.
        for b in 'ACGT':
            if b != cb[1]:
                altcodon = list(codon)
                altcodon[frame] = b
                altcodon = ''.join(altcodon)
                a_aa = translate[altcodon]
                if a_aa != refaa:
                    #its a nonsynonymous change.
                    syndata[cb]['non'].append(b)
                else:
                    syndata[cb]['syn'].append(b)
    return syndata
    
good = 0
bad = 0
sitesyn = {}
#the above will have chrocoord key and value of the synonymity of the mutations at the site.
for mrna in db.features_of_type("mRNA", order_by='start'):
    full_seq = []
    for cds in db.children(mrna, featuretype='CDS',order_by='start'):
        seq = cds.sequence('riv-pilon.fa', use_strand = True) #reverse complemented if strand is negative
        #assign each location in the sequence a coordinate, then remove ones from the start/end based on phase
        scv = [(c,b) for c,b in zip(list(range(cds.start, cds.end+1)),seq)] #scv stands for sequence coordinate vector
        #crop from the start/end based on phase.
        if cds.strand == '+':
            scv = scv[int(cds.frame):]
        elif cds.strand == '-':
            scv = scv[:-int(cds.frame)]
        #else:
            #print("We Got Problems.")
        full_seq.extend(scv)
    #now we proceed to get our synonymity stuff
    #lets track how many of these are actually the correct length
    #I bet its 1/3, because there's no fucking way you'd expect a stitched CDS to actually be composed of fucking CODONS
    try:
        assert len(full_seq)%3 == 0
        good += 1
    except:
        #print("Combined CDS sequence is the wrong length. of course.")
        #print(mrna)
        #print(len(full_seq), len(full_seq)%3)
        bad += 1
        continue
    sd = get_synonymity_info(full_seq)
    #sitesyn.update({(cds.chrom, k):v for k,v in sd.items()}) #add chromosome to key.
    #for each site in sd, print the data out in a format I can parse to split the pileup up later.
    for k,v in sd.items():
        if k[1] != 'N':
            data = [str(v) for v in [cds.chrom, k[0], k[1], ','.join(v['syn']), ','.join(v['non'])]]
            print('\t'.join(data))