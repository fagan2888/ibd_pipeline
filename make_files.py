#!/bin/usr/env python3
from datetime import datetime
import numpy as np
import os
import sys
from input_code import LoadFiles
from parse_code import FileParser

class FileConvert(object):
    """FileConvert() accepts the files that have been validated. The output will be map files with genetic position and plink phased ped files
    """

    def start_pedfile(sample_file):
        """start_pedfile() returns the first five columns of the pedfile
        will need to modify when figuring out how to accomodate for variation in sample files
        """
        pedfile_start=sample_file[:, [0,1,3,4,5,6]]
        return(pedfile_start)

    
    def make_gpos_mapfile(haps_file, hapmap, outfile):
        """make_gpos_mapfile returns (something) that can be saved as a mapfile with genetic position
        output style: [chrom] [rs] [genetic pos] [pos]
        """

        mappos = list()
        mapgpos = list()

        chromin = haps_file[:, 0].tolist()
        rsin = haps_file[:, 1].tolist()
        posin = haps_file[:, 2].tolist()

        line = hapmap.readline()
        line = hapmap.readline()
        while line:
            line = line.strip().split()
            pos = int(line[1])
            gpos = float(line[3])
            mappos.append(pos)
            mapgpos.append(gpos)
            line = hapmap.readline()

        index1 = 0
        index2 = 0
        while index1 < len(posin):
            pos = int(posin[index1])
            rs = rsin[index1]
            chrom = chromin[index1]
            if pos == mappos[index2]:
                ##the 1000 Genomes site was genotopyes as part of the map (comment directly from original code)
                outfile.write(' '.join([chrom, rs, str(mapgpos[index2]), str(pos)]) + '\n')
                index1 = index1 + 1
            elif pos < mappos[index2]:
                ## current position in interpolation before marker
                if index2 == 0:
                    ## before the first site in the map (genetic position = 0)
                    outfile.write(' '.join([chrom, rs, str(mapgpos[index2]), str(pos)]) + '\n')
                    index1 = index1 + 1
                else:
                    ## interpolate
                    prevg = mapgpos[index2 - 1]
                    prevpos = mappos[index2]
                    frac = (float(pos) - float(mappos[index2 - 1]))/ (float(mappos[index2]) - float(mappos[index2 - 1]))
                    tmpg = prevg + frac* (mapgpos[index2] - prevg)
                    outfile.write(' '.join([chrom, rs, str(tmpg), str(pos)]) + '\n')
                    index1 = index1 + 1
            elif pos > mappos[index2]:
                ## current position in iterpolation after marker
                if index2 == len(mappos) - 1:
                    ## after the last site in the map (genetic position = maximum in map, note could try to extrapolate based on rate instead)
                    outfile.write(' '.join([chrom, rs, str(mapgpos[index2]), str(pos)]) + '\n')
                    index1 = index1 + 1
                else:
                    ## increment the marker
                    index2 = index2 + 1
                                                                      
    
    def convert_haps(haps_file):
        """convert_haps() takes the haps_file and converts the binary section (0,1) to (A,T,C, or G)
        I want to optimize this there has got to be a better way
        """
        haps_list = haps_file.tolist()
        hapslist = []
        for line in haps_list:
            base1  = line[3]
            base2 = line[4]
            nums = line[5:]
            nums = map(lambda x: x.replace('0', base1), nums)
            nums = map(lambda x: x.replace('1', base2), nums)
            hapslist.append(list(nums))
        return(hapslist)


    def make_pedfile(hapslist, pedlen, pedfile_start):
        lenhaps = int(len(hapslist))
        t_hapslist = np.array(hapslist)
        r_hapslist = t_hapslist.reshape(lenhaps,pedlen,2)
        c_hapslist = np.concatenate(r_hapslist[0:], axis=1)
        final_list = np.append(pedfile_start, c_hapslist, axis=1)
        return(final_list)


## Load files and uses class from input_code
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
FileParser.validate_binary(haps_file)
stop_ps = datetime.now()
time_diff_ps = stop_ps - start_ps
print("Parsing your files takes {}".format(time_diff_ps))


## File convert section
start_fc = datetime.now()
pedfile_start = FileConvert.start_pedfile(sample_file)   
pedlen = int(len(pedfile_start))

hapmap = open(sys.argv[3])

outfile = open(sys.argv[4], "w")
FileConvert.make_gpos_mapfile(haps_file, hapmap, outfile) 

hapslist = FileConvert.convert_haps(haps_file)

final_list = FileConvert.make_pedfile(hapslist, pedlen, pedfile_start)

np.savetxt(sys.argv[5], final_list, fmt='%s', delimiter=' ')
stop_fc = datetime.now()
time_diff_fc = stop_fc - start_fc
print("Making the map file and plink ped file takes {}".format(time_diff_fc))
