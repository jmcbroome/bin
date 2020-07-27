def sdconvert(smat_f, out):
    '''
    Function that converts a sparse matrix (as a list of string entries) to a dense matrix.

    This is pretty slow but it really do be like that sometimes.
    '''
    import numpy as np
    smat = []
    with open(smat_f, 'r') as sparse:
        maxn = 0
        for entry in sparse.readlines():
            splent = [int(e) for e in entry.strip().split()] 
            if splent != []:
                smat.append(splent)
                if splent[0] > maxn:
                    maxn = splent[0]
    dmat = np.ones(shape = (maxn+1,maxn+1))
    for ent in smat:
        dmat[ent[0]][ent[1]] = int(ent[2]) + 1
        dmat[ent[1]][ent[0]] = int(ent[2]) + 1
    with open(out, 'w+') as dense:
        for rind, row in enumerate(dmat):
            for cind, col in enumerate(row):
                print(dmat[rind][cind], end = '\t', file = dense)
            print('', end = '\n', file = dense)

def main():
    # smatf = input()
    # out = input()
    smatf = 'NWM9_1.dat'
    out = smatf[:-4] + '_dense.txt'
    sdconvert(smatf, out)

main()