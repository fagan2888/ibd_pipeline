#!/bin/usr/env python3
from datetime import datetime
import numpy as np
import os
import sys
from input_code import LoadFiles

class FileParser(object):
    """FileParser() class parses(is this right?) the sample file and haplotype file.
    """

    def check_same_sample(sample_file,haps_file):
        """check_same_sample() tests for consistancy between the files (2 haplotypes per person)
        """
        sample_shape = sample_file.shape
        haps_shape = haps_file.shape
        if sample_shape[0] == (haps_shape[1] - 5)/2:
            print("Each person has 2 haplotypes")
        else:
            raise ValueError("Each person does not have 2 haplotypes.") ## Is this the right kind of error?


    def check_alleles(haps_file):
        """Ensures that the reference and alternative alleles are different.
        """
        diff_alleles = np.array(haps_file[:,3] == haps_file[:,4])
        if True in diff_alleles:
            raise ValueError("Your reference and alternative alleles should not be the same.")
        else:
            print("All of your reference and alternative alleles are not the same.")


    def nucleotide_test(nucleotide):
        """Tests if the nucleotides are A,C,T, or G.
        """
        if nucleotide.upper() in ['A', 'C', 'T', 'G']:
            return(True)
        else:
            return(False)
    
    def validate_nucleotides(array):
        """Validates results from nucleotide_test.
        """
        if False in array:
            raise ValueError("The reference or alternative alleles should be A, T, C, or G") ## Ask how to customize this ( I want to say reference or alternative and include loc or false)
        else:
            print("All reference or alternatives alleles are correct.")

    def validate_binary(haps_file):
        """Ensures that the column 6 to the end of the file is 0 or 1.
        """
        log_array_haps = np.logical_or(haps_file[:,5:-1] == '0', haps_file[:,5:-1] =='1')
        if False in log_array_haps:
            raise ValueError("This is not a properly formated Haplotype file. The haplotypes should be 0's and 1's.")
        else:
            print("Binary representations of SNPS are correct.")


## Loads files and uses class from input_code
start_lf = datetime.now()
sample_file = LoadFiles.load_sample(sys.argv[1])
haps_file = LoadFiles.load_haps(sys.argv[2])
stop_lf = datetime.now()
time_diff_lf = stop_lf - start_lf
print("Loading your files takes {}".format(time_diff_lf))  

## Parsing section
start_ps = datetime.now()
FileParser.check_same_sample(sample_file, haps_file)
FileParser.check_alleles(haps_file)

nucfunc = np.vectorize(FileParser.nucleotide_test)
log_array_ref = nucfunc(haps_file[:,3].astype(str))
log_array_alt = nucfunc(haps_file[:,4].astype(str))

FileParser.validate_nucleotides(log_array_ref)
FileParser.validate_nucleotides(log_array_alt)
stop_ps = datetime.now()
time_diff_ps = stop_ps - start_ps
print("Parsing your files takes {}".format(time_diff_ps))




