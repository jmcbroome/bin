import re
def motif_finder(f1, motif, f2 = 'motifs.txt'):
    #first, you have to write the column names to the file
    with open(f2, 'w+') as file_out:
        column_string = #fill this in- remember \t is the symbol for a tab space
        print(column_string, file = file_out)
        #now read the data file, and iterate over all rows in it
        with open(f1, 'r') as file_in:
            for i, line in enumerate(file_in): 
                #enumerate is a helpful function; basically it takes a list [a,b,c] and turns it into [(1,a),(2,b),(3,c)] which makes tracking indexes easier. then I can assign both the 1 and the a to two variables with the for i, entry syntax.
                if line[0] != '>':
                    #for whatever reason, it looks like your teacher wants you to replace the name "Protein1" with "Seq1"?
                    name = #fill this in too! you're probably going to want to use that i I defined for this. Or you can write something else, whatever works!
                    #they also want you to count hits, which is the number of instances of the motif you find, presumably.
                    hits = re.findall(motif, line) #re.findall() actually returns a list, which has 1 entry for each place where it finds the motif you give it
                    hit_number = len(hits)
                    #there are two different ways to handle prepping the actual entry for printing, you can either write it out with a \t between every name and combine, or you can put them all in a list and use '\t'.join(thelist)
                    entry = 
                    print(entry, file = file_out) #default will include a newline \n so they won't be all on the same line

                    print('\t'.join(['Seq'+str(i),motif,len(re.findall(motif,line)),','.join(re.findall(motif,line))]),file = file_out)