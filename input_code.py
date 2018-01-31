#!/bin/usr/env python3
from datetime import datetime
import numpy as np
import os
import sys


class LoadFiles(object):
    """LoadFiles
    Class for loading in file. Version (##) will load in sample files and haplotype files.
    Version (## + 1) will load in plink format pedigree files and map files.
    Both versions require the files to be phased. 
      

    """
    def load_sample(path):
        """load_sample() Takes in the sample file as input and this file has two rows that make
        up the header. They will be skipped.
        Sample File Format:
        1: ID_1 
        2: ID_2
        3: missing
        4-X: covariates
        
        """
        if os.path.isfile(path):
            sample_file = np.loadtxt(path, dtype=str, delimiter=' ', skiprows=2)
            # FR: will this use the first line as the header?
            return(sample_file)
        else:
            raise IOError("Sample File does not exist.")

    def load_haps(path):
        """load_haps() Takes in the haplotype file as input and this file format should be
        1: Chromosome #
        2: RSID
        3: Position
        4: Reference Allele
        5: Alternative Allele
        6-X: Binary representations of the alleles 
        Note: 0 maps to the reference allele and 1 maps to the alternative allele
        """
        if os.path.isfile(path):
            haps_file = np.loadtxt(path, dtype=str, delimiter=' ')
            return(haps_file)
        else:
            raise IOError("Haplotype File does not exist.")

"""
start = datetime.now()

sample_file = LoadFiles.load_sample(sys.argv[1])
haps_file = LoadFiles.load_haps(sys.argv[2])

stop = datetime.now()
time_diff = stop - start
print(time_diff)
"""
