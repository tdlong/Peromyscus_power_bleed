import sys
import numpy as np

first_iteration = True
for line in sys.stdin:
    #initialize line, removing CHROM and POS columns
    dosages = np.delete(line.split(), [0, 1]).astype("float32")
    
    # create Gb on first iteration
    if first_iteration:
        Gb = np.zeros(len(dosages), dtype = "float32")
        first_iteration = False
        
    #add current dosages to Gb (which is a sum of dosages)
    np.add(Gb, dosages, out = Gb)
    
######################

individuals = np.array(open("/share/adl/pnlong/mouseproject/save_data/all_individuals.tsv", "r").read().split())

scale = lambda vector: np.array([((x - np.mean(vector)) / np.std(vector)) for x in vector], dtype = "float32")
Gb = np.sqrt(45) * scale(Gb)

#print it out
print("INDIVIDUAL", "GBV_VALUE", sep = "\t", end = "\n")
for i in range(len(individuals)):
    print(str(individuals[i]), str(Gb[i]), sep = "\t", end = "\n")