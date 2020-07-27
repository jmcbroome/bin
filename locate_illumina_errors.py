#the goal of this script is to parse in a sam file and return a set of locations and sequence information for likely illumina errors
#these are bases which are N within a consensus while bases to either side of them are non-N.
import sys

#takes a sam file of consensus sequences
#prints a table output that is read name, chromosome, coordinate, prior, latter, tab separated
#identity of the error'd base isn't derivable from the sam, would need to be incorporated into the pileup step, going to pass on that for now
#though you could probably look it up in the reference and it'd probably be that.

print('\t'.join(['ReadName','Chro','Loc','Prior','Latter','ConQual', 'Motif']))
for entry in sys.stdin:
    read_name, _, chro, loc, _, _, _, _, _, seq, quals = entry.strip().split()
    #iterate through the sequence string looking for solo N's.
    for i, bq in enumerate(zip(seq, quals)):
        b,q = bq
        if b == 'N' and q != '0':
            #an illumina error! bonus encoding- S = start of sequence, E = end of sequence.
            if i == 0:
                prior = 'S'
            else:
                prior = seq[i-1]
            if i == len(seq)-1:
                latter = 'E'
            else:
                latter = seq[i+1]
            if prior == "N" or latter == 'N':
                continue
            #new column- collect the 3mer preceding motif, which is apparently important to illumina errors. See if I can replicate enrichment of CGG and GGG motifs.
            motif = seq[i-3:i]
            if any([m == 'N' for m in motif]):
                continue #once again, skip these.
            print('\t'.join([read_name, chro, str(int(loc) + i), prior, latter, q, motif]))
