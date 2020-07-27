#!/usr/bin/env python3

import sys

#set id conversion dictionary.
# convert = {
#     "NC_004354.4":"chrX",
#     "NT_033779.5":"chr2L",
#     "NT_033778.4":"chr2R",
#     "NT_037436.4":"chr3L",
#     "NC_004252.4":"chr4",
#     "NT_033777.3":"chr3R",
#     "NC_024512.1":"chrY"
# } #Fly chromosomes.

convert = {
    'NC_000001.11':'chr1',
    'NC_000002.12':'chr2',
    'NC_000003.12':'chr3',
    'NC_000004.12':'chr4',
    'NC_000005.10':'chr5',
    'NC_000006.12':'chr6',
    'NC_000007.14':'chr7',
    'NC_000008.11':'chr8',
    'NC_000009.12':'chr9',
    'NC_0000010.11':'chr10',
    'NC_0000011.10':'chr11',
    'NC_0000012.12':'chr12',
    'NC_0000013.11':'chr13',
    'NC_0000014.9':'chr14',
    'NC_0000015.10':'chr15',
    'NC_0000016.10':'chr16',
    'NC_0000017.11':'chr17',
    'NC_0000018.10':'chr18',
    'NC_0000019.10':'chr19',
    'NC_0000020.11':'chr20',
    'NC_0000021.9':'chr21',
    'NC_0000022.11':'chr22',
    'NC_0000023.11':'chrX',
    'NC_0000024.10':'chrY',
} #human chromosomes.

for entry in sys.stdin:
    spent = entry.strip().split()
    nname = convert.get(spent[0], None)
    if nname != None:
        print('\t'.join([nname] + spent[1:]))