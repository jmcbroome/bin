import sys
import math
for entry in sys.stdin:
    spent = entry.strip().split()
    #divide it into 5k chunks, rounding off the edges.
    #first need to ceil and floor the boundary edges to the nearest 5k
    cstart = math.ceil(int(spent[1])/5000)*5000
    cend = math.floor(int(spent[2])/5000)*5000
    #create a range
    chunks = list(range(cstart,cend,5000))
    for i in range(1,len(chunks)):
        #assemble a new entry
        nent = [spent[0], str(chunks[i-1]), str(chunks[i])] + spent[3:]
        print('\t'.join(nent)) #print to stdout.
