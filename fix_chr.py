#!/usr/bin/env python3

import sys
change = {'0':"2L", '1':"2R", '2':"3L", '3':"3R", '4':'4', '5':"M", '6':'X'}
for entry in sys.stdin:
    if entry[0:2] == "##":
        spent = entry.strip().split(',')
        if spent[0][-1] in change:
            spent[0] = spent[0][:-1] + change[spent[0][-1]]
        print(','.join(spent))
    else:   
        spent = entry.strip().split()
        if spent[0] in change:
            spent[0] = change[spent[0]]
        print('\t'.join(spent))