# README
# This program averages the Mjj Values for the previously calculated sum values for the kinship matrix
# Phillip Long
# February 20, 2020

import sys
import numpy as np

individuals = tuple(open("/share/adl/pnlong/mouseproject/save_data/all_individuals.tsv", "r").read().split())

## determine crosses making a 81010 x 2 array
crosses = []
for individual1 in range(0, len(individuals) - 1): #so that when it gets to the second last row (last iteration), individual2 will be the last row
    for individual2 in range((individual1 + 1), len(individuals)):
        cross = (individuals[individual1], individuals[individual2])
        crosses.append(cross)
##
del individuals

# print the header
print("INDIVIDUAL.x", "INDIVIDUAL.y", "MJJ", sep = "\t")

counter = 0
# MAIN READER
# for loops is being used because to calculate the Mjj_sum,
# in an ideal world, the counter and Mjj_sum would come out as one line
# but in reality, that would take up too much memory,
# so we had to split up the task into an array job;
# thus, the are many different Mjj_sum and counter values that
# need to be summer even further.
for line in sys.stdin:
    row = np.array(line.split()) # create the row as numpy array
    
    if counter == 0: # create Mjj_sum on first iteration
        Mjj_sum = np.zeros((len(row) - 1), dtype = "float32") # (len(row) - 1) because first index is the counter value
        
    read_count = int(row[0])
    row = np.delete(row, 0).astype("float32")
    
    counter += read_count
    
    np.add(Mjj_sum, row, out = Mjj_sum)

###############################
# END OF CODE TO CALCULATE SUMS
Mjj_avg = np.array((Mjj_sum / float(counter)), dtype = "float32")

# END OF CODE TO CALCULATE Averages
###############################

#print out the info
for i in range(len(crosses)):
    print(*crosses[i], Mjj_avg[i], sep = "\t", end = "\n")