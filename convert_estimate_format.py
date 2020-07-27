#!/usr/bin/env python3

def convert():
    chromlen = {'chr2L':23515712, 'chr2R':25288936, 'chr3L':25288936, 'chr3R':32081331, 'chrX':23544271} #values as of DM6 assembly, chr 2/3/X only for first testing
    with open('correctBreakpointEstimates.txt', 'r') as original:
        with open("bpest.txt", 'w+') as new:
            for line in original.readlines():
                nline = line.strip().split()
                if nline[1] == 'X':
                    print('chrX\t' + str(nline[3]) + '\t' + str(nline[4]), file = new)
                if len(nline) == 6:
                    if nline[1][1:] == 'L':
                        print('chr' + str(nline[1]) + '\t' + str(nline[3]) + '\t' + str(nline[4]), file = new)
                    elif nline[1][1:] == 'LR':
                        print('chr' + nline[1][0] + 'L\t' + str(nline[3]), file = new)
                        if int(nline[4]) - chromlen['chr' + nline[1][0] + 'L'] > 0:
                            ncoord = int(nline[4]) - chromlen['chr' + nline[1][0] + 'L']
                            print('chr' + nline[1][0] + 'R\t' + str(ncoord), file = new)
                        else:
                            print('chr' + nline[1][0] + 'R\t' + nline[4], file = new)
                    elif nline[1][1:] == 'R':
                        ncoord1 = int(nline[3]) - chromlen['chr' + nline[1][0] + 'L']
                        if ncoord1 < 0:
                            ncoord1 = int(nline[3])
                        ncoord2 = int(nline[4]) - chromlen['chr' + nline[1][0] + 'L']
                        if ncoord2 < 0:
                            ncoord2 = int(nline[4])
                        print('chr' + nline[1][0] + 'R\t' + str(ncoord1) + '\t' + str(ncoord2), file = new)
                else:
                    if nline[1][1:] == 'L':
                        print('chr' + str(nline[1]) + '\t' + str(nline[2]) + '\t' + str(nline[3]), file = new)
                    elif nline[1][1:] == 'LR':
                        print('chr' + nline[1][0] + 'L\t' + str(nline[2]), file = new)
                        if int(nline[3]) - chromlen['chr' + nline[1][0] + 'L'] > 0:
                            ncoord = int(nline[3]) - chromlen['chr' + nline[1][0] + 'L']
                            print('chr' + nline[1][0] + 'R\t' + str(ncoord), file = new)
                        else:
                            print('chr' + nline[1][0] + 'R\t' + nline[3], file = new)
                    elif nline[1][1:] == 'R':
                        ncoord1 = int(nline[2]) - chromlen['chr' + nline[1][0] + 'L']
                        ncoord2 = int(nline[3]) - chromlen['chr' + nline[1][0] + 'L']
                        if ncoord1 < 0:
                            ncoord1 = int(nline[2])
                        if ncoord2 < 0:
                            ncoord2 = int(nline[3])
                        print('chr' + nline[1][0] + 'R\t' + str(ncoord1) + '\t' + str(ncoord2), file = new)
convert()