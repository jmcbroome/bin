#!/usr/bin/env python3

#import
import argparse

#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    # parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-b', '--bed', help = 'path to bed file to correct score values of. Score is final column.')
    parser.add_argument('-c', '--context', help = 'path to alternate bed file with equivalent intervals and the score to subtract')
    parser.add_argument('-o', '--out', help = 'name of corrected file.')
    args = parser.parse_args()
    return args

def read_cluster(path):
    #read the cluster score data into a dictionary with ordered array values.
    clusterd = {}
    with open(path) as inf:
        for entry in inf:
            chro, start, end, size, score = entry.strip().split()
            size = int(size)
            start = int(start)
            end = int(end)
            try:
                score = float(score)
            except ValueError:
                continue
            if size not in clusterd:
                clusterd[size] = {}
            if chro not in clusterd[size]:
                clusterd[size][chro] = []
            clusterd[size][chro].append((start, end, score))
    #order the lists for each cluster size for each chromosome.
    for size,chros in clusterd.items():
        for c, locs in chros.items():
            clusterd[size][c] = sorted(locs, key = lambda x:x[0])
    return clusterd
  
def make_locd(tuples):
    locd = {}
    for start, end, score in tuples:
        locd[start] = (end,score)
    return locd

def read_score(path):
    scored = {}
    with open(path) as inf:
        for entry in inf:
            chro, start, end, score = entry.strip().split()
            start = int(start)
            end = int(end)
            score = float(score)
            if chro not in scored:
                scored[chro] = []
            scored[chro].append((start, end, score))
    return scored

def compare_bed(bed1, bed2):
    nbed = {}
    for size, chros in bed1.items():
        if size not in nbed:
            nbed[size] = {}
        for c, locs in chros.items():
            if c not in nbed[size]:
                nbed[size][c] = []
            b1d = make_locd(locs)
            try:
                b2d = make_locd(bed2[size][c])
            except KeyError:
                print("Chromosome {} missing context score data".format(c))
                continue
            for start, endscore in b1d.items():
                try:
                    nscore = endscore[1] - b2d[start][1]
                except KeyError:
                    continue
                nbed[size][c].append((start, endscore[0], nscore))
    return nbed

def write_bed(clusterd, outf):
    with open(outf, 'w+') as out:
        for size, chros in sorted(clusterd.items()):
            for c, locs in sorted(chros.items()):
                for l in sorted(locs):
                    start, end, score = l
                    line = '\t'.join([str(c),str(start),str(end),str(size),str(score)])
                    print(line,file = out)

def main():
    args = argparser()
    clusters = read_cluster(args.bed)
    contexts = read_cluster(args.context)
    corrected = compare_bed(clusters,contexts)
    write_bed(corrected, outf = args.out)

if __name__ == "__main__":
    main()