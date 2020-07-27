#!/usr/bin/python
from scipy.optimize import minimize 
from sys import argv
from math import log10 
import csv
import gzip
import random

### initial prediction
breakpoint = [ float( argv[2] ), float( argv[3] ) ]

### parameters
distance = 5e6 					#### distance to the predicted breakpoints for read pairs we consider
proportion_retained = float(argv[4])

### data objects
p1 = []
p2 = []

### compute distance function
def compute_distance( breakpoint ) : 

	### take in teh values
	pos1, pos2 = breakpoint

	### store distance here 
	dist = 0 

	### iterate through positions
	for i in range( len(p1)-1 ) :

		## now go through and update our lines and compute the distance for the points
		if ( ( p1[i] < pos1 and p2[i] > pos2 ) or ( pos1 < p1[i] < pos2 and pos1 < p2[i] < pos2 ) or ( p1[i] <= p2[i] < pos1 ) ): ## ( p1[i] > pos2 and p2[i] > pos2 ) ) : 
			dist += ( p2[i] - p1[i] ) 
		if ( p1[i] < pos1 and pos1 < p2[i] < pos2 ) :
			dist += ( pos1 + ( pos2 - p2[i] ) - p1[i] ) 
		elif ( pos1 < p1[i] < pos2 and pos2 < p2[i] ) : 
			dist += ( p2[i] - ( pos2 - ( p1[i] - pos1 ) ) )

	### return the total distance spanned by the reads 
	return dist 

### read data into paired-lists 
with gzip.open( argv[1] ) as tsv :

	### split on tab
    for line in csv.reader(tsv, delimiter="\t") :

    	### check to make sure we're close enough to consider the read
    	if ( ( abs( float(line[1]) - breakpoint[0] ) < distance or abs( float(line[1]) - breakpoint[1] ) < distance ) and ( abs( float(line[2]) - breakpoint[0] ) < distance or abs( float(line[2]) - breakpoint[1] ) < distance ) ):

    	### alternatively, only consider sites where each read is near either breakpoint 
    	##if ( abs( float(line[1]) - breakpoint[0] ) < distance and abs( float(line[2]) - breakpoint[1] ) < distance ):


    		### check uniqueness of read-pairs 
    		if ( p1[-1:] == float(line[1]) and p2[-1:] == float(line[2]) ):
    			next()

		### subsample read pairs
		elif ( random.random() < proportion_retained ) :

	    		### append read positions to list 
    			if ( float( line[1] ) < float( line[2] ) ) :
	    			p1.append( float( line[1] ) )
    				p2.append( float( line[2] ) )

    			### or if alternative, add in other order
       			elif ( float( line[1] ) > float( line[2] ) ) :
	    			p1.append( float ( line[2] ) )
    				p2.append( float ( line[1] ) )  			

### now do the optimization
estimate = minimize(compute_distance, breakpoint, method="Nelder-Mead", options={'maxiter':5000,'maxfev':5000} )
print estimate.fun, estimate.nfev, estimate.success, estimate.x

#for i in range ( 1, 100 ) :
#	for j in range ( 1, 100 ) : 
#		print breakpoint[0]+i*1000, breakpoint[1]+j*1000, compute_distance( [breakpoint[0]+i*1000,breakpoint[1]+j*1000] ) 


