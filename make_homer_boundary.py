#process a 2D bed file of single boundary lines into regions of 1k bases around each line.
import sys
for entry in sys.stdin:
    chro, start, stop, null, score = entry.strip().split()
    print('\t'.join([chro, str(int(start)-500), str(int(start)+500)]))
    print('\t'.join([chro, str(int(stop)-500), str(int(stop)+500)]))
