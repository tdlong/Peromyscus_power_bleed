# README
# This program calculates the unaveraged Mjj values (the sums) for the kinship matrix
# Phillip Long
# February 20, 2020

# python Mjj_sum.py "/share/adl/pnlong/mouseproject/save_data/individuals.tsv" "/share/adl/pnlong/mouseproject/save_data/all_individuals.tsv"
# sys.argv[1] = path to subsetted list of individuals
# sys.argv[2] = path to full list of individuals

import sys
import numpy as np
from math import exp


# import lists of all individuals and individuals needed to extract
individual_IDs_full = tuple(open(str(sys.argv[2]), "r").read().split())
individual_IDs = tuple(open(str(sys.argv[1]), "r").read().split())
N = len(individual_IDs)
number_of_crosses = int(N * (N - 1) * 0.5)


# get the indices of the individuals that we need to extract
original_positions = np.array(tuple(individual_IDs_full.index(ID) for ID in individual_IDs))
starting_indices_to_extract = tuple(map(lambda pos: 8 * pos, original_positions))
indices_to_extract = np.array(list(map(lambda pos: tuple(range(pos, pos + 8)), starting_indices_to_extract)), dtype = np.int).flatten()


Mjj_sum = np.zeros(shape = number_of_crosses, dtype = "float32")


calculate_Mjj = lambda dose1, dose2: (dose1[0] * dose2[0]) + (0.5 * ((dose1[0] * dose2[1]) + (dose1[1] * dose2[0]) + (dose1[1] * dose2[1])))
# (1.0 * probability_matrix[row1][dose][0] * probability_matrix[row2][dose][0]) + (0.5 * probability_matrix[row1][dose][0] * probability_matrix[row2][dose][1]) + (0.5 * probability_matrix[row1][dose][1] * probability_matrix[row2][dose][0]) + (0.5 * probability_matrix[row1][dose][1] * probability_matrix[row2][dose][1])
    
    
counter = 0
#main reader
for line in sys.stdin:
    counter += 1
    
    # split and create the row as a numpy array,
    # np.delete the first two columns (0, 1) --> (CHROM & POS),
    # then cast as float32
    row = np.delete(np.array(line.strip().split()), (0,1)).astype("float32")
    row = np.array(tuple(row[index] for index in indices_to_extract)) # filter down to just the individuals we want
    
    row2D = row.reshape(N, 8) # set data type of row to precise float and turn into matrix: number of individuals by number of haplotypes
    # new 2D list for of the row, delete "row" variable to save memory
    del row
    
    probability_matrix = np.ndarray(shape = (N, 8, 2), dtype = "float32")
    
    # calculate 2 probabilities for each individual and {8} haplotype combination
    for row in range(N): #len(row2D) = N
        for col in range(8): #len(row2D[row]) = 8
            
            # calculate the denominator
            denominator = exp(-(row2D[row][col] - 0)**2 / 0.18) + exp(-(row2D[row][col] - 1)**2 / 0.18) + exp(-(row2D[row][col] - 2)**2 / 0.18)
            
            # calculate two values for each element of row2D
            probability_matrix[row][col][0] = exp(-(row2D[row][col] - 2)**2 / 0.18) / denominator
            probability_matrix[row][col][1] = exp(-(row2D[row][col] - 1)**2 / 0.18) / denominator
            
    del row2D # it is not needed anymore; to save memory
    
    
    
    Mjj_vector = np.zeros(shape = number_of_crosses, dtype = "float32")
    
    i = 0
    for row1 in range(0, (len(probability_matrix) - 1)):
        for row2 in range((row1 + 1), (len(probability_matrix))):
            Mjj_vector[i] = sum(map(calculate_Mjj, probability_matrix[row1], probability_matrix[row2]))
            i += 1

    del probability_matrix # save memory
 
        
    #add current Mjj values to the sum
    np.add(Mjj_sum, Mjj_vector, out = Mjj_sum)

###############################
# print it out
# print the counter to help in averaging in the next step
print(counter, *Mjj_sum, sep = "\t")