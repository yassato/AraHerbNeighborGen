####################################
# Subsetting the full imputed SNPs #
####################################

# load libraries
import h5py
import pandas as pd
import numpy as np
import csv

# load data
# The full_imputed_SNP is available at https://aragwas.1001genomes.org/#/download-center
# 1. Access via "For other genotype versions, please refer to this google drive page" (accessed on 26 Oct 2022).
# 2. A compressed file "full_imputed_SNP_MATRIX_2Jun2022.tar.gz" includes "all_chromosomes_binary.hdf5" used below.
# 3. As recommended, the input file .hdf5 was handled using Python below
input_file = "../data/all_chromosomes_binary.hdf5"
h5file = h5py.File(input_file,"r")

snps = h5file["snps"].value
accs = h5file["accessions"].value
pos = h5file["positions"].value

acID = []
for i in accs:
    acID.append(int(i))

acID = np.array(acID)
gwasID = pd.read_csv("../data/gwasIDlist.csv") # load gwasIDs

acc_list = []
for i in gwasID["GWASid"]:
    place = np.where(acID == i) # search a focal accession
    acc_list.append(int(place[0]))
    
# export the subset data
sub_snps = snps[:,acc_list]
del snps
sub_snps = pd.DataFrame(sub_snps)
sub_snps.to_csv("../data/sub_snps.csv") # export SNPs

pos = pd.DataFrame(pos)
pos.to_csv("../data/positions.csv") # export positions
