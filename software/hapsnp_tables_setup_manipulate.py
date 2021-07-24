# README
# this program wrangles freshly extracted vcf files into a more tidy .tsv format
# Phillip Long
# December 17, 2020

# python /share/adl/pnlong/mouseproject/software/hapsnp_tables_setup_manipulate.py "snp" "/share/adl/pnlong/mouseproject/save_data/individuals.tsv" "/share/adl/pnlong/mouseproject/save_data/all_individuals.tsv"
# python /share/adl/pnlong/mouseproject/software/hapsnp_tables_setup_manipulate.py "hap" "/share/adl/pnlong/mouseproject/save_data/individuals.tsv" "/share/adl/pnlong/mouseproject/save_data/all_individuals.tsv"
# sys.argv[1] = statistical test type of table (either "snp" or "hap")
# sys.argv[2] = path to subsetted list of individuals
# sys.argv[3] = path to full list of individuals


import sys
import numpy as np

# import lists of all individuals and individuals needed to extract
position_column = "POS"
individual_IDs_full = tuple([position_column, ] + open(str(sys.argv[3]), "r").read().split())
individual_IDs = tuple(open(str(sys.argv[2]), "r").read().split())

# print the header, get indices needed to extract
table_type = str(sys.argv[1])
if table_type == "snp":
    print("POS", *individual_IDs, sep = "\t", end = "\n")
    # get the indices of the individuals that we need to extract
    indices_to_extract = tuple(individual_IDs_full.index(ID) for ID in (position_column, *individual_IDs))
elif table_type == "hap":
    column_list = np.array(tuple(map(lambda ID: tuple(f"HAP_{i}-{ID}" for i in range(1, 9)), individual_IDs))).flatten()
    print("POS", *column_list, sep="\t", end="\n")
    
    # get the indices of the individuals that we need to extract
    original_positions = np.array(tuple(individual_IDs_full.index(ID) for ID in individual_IDs))
    
    # 8 because 8 values per haplotype
    # (x - 1) because if the input is 32, it is actually the 31st group of 8 (due to the position column originally causing a +1)
    # +1 at the end because it should be +2... due to
    # - +1 from position column and
    # - +1 from the fact that the number you are getting from the multiply by 8 is the 8th value of the previous haplotype, not the 1st value of the current haplotype
    # but -1 from +2 because python starts counting at 0
    # so if you are starting your 2nd haplotype (normally should be the 9th index),
    # it would return 9 because even though it SHOULD be 1 less because of python starting at 0, the position column at the start makes up for that
    transform_position = lambda x: (8 * (x - 1)) + 1
    starting_indices_to_extract = tuple(map(transform_position, original_positions))
    
    indices_to_extract = np.insert(np.array(list(map(lambda x: tuple(range(x, x + 8)), starting_indices_to_extract)), dtype = np.int).flatten(), 0, 0)
else:
    print("Faulty input for the table type.")

# print out bulk of the file
i = 0
for line in sys.stdin:
    i += 1
    if table_type == "hap" and i % 10 != 1: # if a haplotype_table, only print every 10th line
        continue
    
    full_row = line.split()
    row = (full_row[index] for index in indices_to_extract)
    print(*row, sep = "\t", end = "\n")