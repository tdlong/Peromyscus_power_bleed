# README
# this program wrangles freshly extracted vcf files into a more tidy .tsv format
# Phillip Long
# July 12, 2021

import sys
import numpy as np

# import lists of all individuals and individuals needed to extract
chromosome_column = "CHROM"
position_column = "POS"
individual_IDs_full = tuple([chromosome_column, position_column, ] + open("/share/adl/pnlong/mouseproject/save_data/all_individuals.tsv", "r").read().split())
individual_IDs = tuple(open("/share/adl/pnlong/mouseproject/save_data/rnaseq_individuals.tsv", "r").read().split())

# get the indices of the individuals that we need to extract
indices_to_extract = tuple(individual_IDs_full.index(ID) for ID in (chromosome_column, position_column, *individual_IDs))

# print out bulk of the file
for line in sys.stdin:
    full_row = line.split()
    row = tuple(full_row[index] for index in indices_to_extract)
    for i in range(len(individual_IDs)):
        # print(CHROM, POS, ID, DNA_GENO)
        print(row[0].replace("Chr", ""), row[1], individual_IDs[i], row[i + 2], sep = "\t", end = "\n")
