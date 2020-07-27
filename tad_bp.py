#5 minute scripting.

def TADRead(filename, binsize = 150000):
    '''
    Reads in a simplified TAD breakpoint matrix (2 column, space/tab delineated, all integers). Argument binsize converts bin coordinates to base coordinates as it reads.
    Returns an array of base values representing outer edges of chromatin domains.
    '''
    tads = []
    with open(filename, 'r') as tadbp:
        for entry in tadbp.readlines():
            ent = entry.strip().split()
            tads.append(ent[0] * binsize)
            tads.append(ent[1] * binsize)
    return set(tads)

def POIRead(filename):
    '''
    Function that translates a file of points of interest into a dictionary format.
    '''
    poidict = {}
    with open(filename) as poi:
        for entry in poi.readlines():
            nent = entry.strip().split()
            poidict[nent[0]] = poidict.get(nent[0], [])  + nent[1:]
    return poidict 

def ttest(tads, poi, chrosize):
    from scipy.stats import ttest_ind
    import random, math
    mins = []
    print("Finding minima.")
    ndist = 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 #effectively infinite, whatever.
    for p in [int(p) for p in poi]:
        for t in [int(t) for t in tads]:
            if p-t > ndist:
                mins.append(ndist)
                break
            else:
                ndist = p-t
        # mins.append(min([abs(int(p)-int(t)) for t in tads]))
    print("Generating distribution.")
    contp = [math.trunc(random.random() * chrosize) for _ in range(100)]
    cmins = []
    print("Finding control minima.")
    ndist = 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 #effectively infinite, whatever.
    for p in contp:
        for t in [int(t) for t in tads]:
            if p-t > ndist:
                cmins.append(ndist)
                break
            else:
                ndist = p-t
    return ttest_ind(mins, cmins, equal_var=False)

def main():
    poifile = 'bpest.txt'
    tadfile = 'results_domains.txt'
    print("Reading POI.")
    poi = POIRead(poifile)
    print("Reading TADs.")
    tads = TADRead(tadfile)
    chromlen = {'chr2L':23515712, 'chr2R':25288936, 'chr3L':25288936, 'chr3R':32081331, 'chrX':23544271} #values as of DM6 assembly, chr 2/3/X only for first testing
    with open("tadnear_out.txt", 'w+') as outf:
        print("Testing.")
        for key, value in chromlen.items():
            print("Chr:{}\tPval:{}".format(key, ttest(tads, poi[key], value)), file = outf)
            print(key + " Finished.")

main()