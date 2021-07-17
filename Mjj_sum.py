# README
# This program calculates the unaveraged Mjj values (the sums) for the kinship matrix
# Phillip Long
# February 20, 2020

import sys
import numpy as np
import math

counter = 0
#main reader
for line in sys.stdin:
    counter += 1
    
    # split and create the row as a numpy array,
    # np.delete the first two columns (0, 1) --> (CHROM & POS),
    # then cast as float32
    row = np.delete(np.array(line.strip().split()), (0,1)).astype("float32")
    
    if counter == 1: # find the number of individuals during the first iteration
        N = len(row) // 8 # N = number of individuals, 8 values per individual's haplotype
    
    row2D = row.reshape(N, 8) # set data type of row to precise float and turn into matrix: number of individuals by number of haplotypes
    # new 2D list for of the row, delete "row" variable to save memory
    del row
    
    probability_matrix = np.ndarray(shape = (len(row2D), len(row2D[0]), 2), dtype = "float32")
    
    # calculate 2 probabilities for each individual and {8} haplotype combination
    for row in range(len(row2D)): #len(row2D) = N
        for col in range(len(row2D[row])): #len(row2D[row]) = 8
            
            # calculate the denominator
            denominator = math.exp(-(row2D[row][col] - 0)**2 / 0.18) + math.exp(-(row2D[row][col] - 1)**2 / 0.18) + math.exp(-(row2D[row][col] - 2)**2 / 0.18)
            
            # calculate two values for each element of row2D
            probability_matrix[row][col][0] = math.exp(-(row2D[row][col] - 2)**2 / 0.18) / denominator
            probability_matrix[row][col][1] = math.exp(-(row2D[row][col] - 1)**2 / 0.18) / denominator
            
    del row2D # it is not needed anymore; to save memory
    
    
    
    Mjj_vector = np.empty(0, dtype = "float32")
    
    for row1 in range(0, (len(probability_matrix) - 1)):
        for row2 in range((row1 + 1), (len(probability_matrix))):
            value = 0.0
            for dose in range(8):
                # sum the combinations of probability values for the 8 haplotype dosages
                value += (1.0 * probability_matrix[row1][dose][0] * probability_matrix[row2][dose][0]) + (0.5 * probability_matrix[row1][dose][0] * probability_matrix[row2][dose][1]) + (0.5 * probability_matrix[row1][dose][1] * probability_matrix[row2][dose][0]) + (0.5 * probability_matrix[row1][dose][1] * probability_matrix[row2][dose][1])
            Mjj_vector = np.append(Mjj_vector, value) # add value after the sum of the 8 values (one per dose in a haplotype) is calculated to Mjj vector
        
    del probability_matrix # save memory
    
    # create the sum array on the first iteration
    if counter == 1:
        Mjj_sum = np.zeros(len(Mjj_vector), dtype = "float32")
        
    #add current Mjj values to the sum
    np.add(Mjj_sum, Mjj_vector, out = Mjj_sum)

###############################
# print it out
# print the counter to help in averaging in the next step
print(counter, *Mjj_sum, sep = "\t")