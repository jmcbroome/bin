import sys
import numpy as np
import math

info = []
last = None
for entry in sys.stdin:
    chro, start, stop, score = entry.strip().split()
    if last != None:
        delta = float(score) - last
        last = float(score)
    else:
        delta = 0
        last = float(score)
    info.append((chro, start, stop, score, delta))
#get the entries which have the lowest deltas (least change from previous)
bgent = sorted(info, key = lambda x:x[4])
for entry in bgent[:math.trunc(len(bgent)/4)]:
    print('\t'.join(entry[:-1]))