#! /usr/bin/env python
#####################################################################
#generate_kmer_distribution.py creates the kmer distribution file needed for est_abundance.py
#Copyright (C) 2016-2019 Jennifer Lu, jlu26@jhmi.edu

#This file is part of Bracken.

#Bracken is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the license, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>.

#####################################################################
#Jennifer Lu, jlu26@jhmi.edu
#01/10/2016
#
#This program evaluates the read distribution of each genome
#and saves the number of reads expected for a given species
#to hit a given taxonomy ID
#
#Note: In this version, error messages have been suppressed. 
#
#Input file: Kraken counts file listing the kraken classification of kmers
#from each genome against the full database. This file must contain the following
#tab-delimited columns:
#   - Identification string for the read
#   - Taxonomy ID of the genome being classified 
#   - Taxonomy ID of the classification kraken assigned to the full read   
#   - Sequence length of the read
#   - Individual Classification Taxonomy IDs for each kmer  
#
#Output file: Kmer distribution file with the following tab-delimited columns:
#   - Classification Taxonomy ID that genomes map to
#   - Genome Taxonomy IDs mapping to that classification (each genome is 
#     space-delimited)
#       - For each genome, the genome taxonomy ID, number of mapped kmers,
#         and number of total kmers are listed, separated by colons ':'
#Methods:
#   - main
#   - parse_single_genome   

import sys, argparse
from time import gmtime
from time import strftime

#parse_single_genome method 
#usage: parses a single line from the input file and extracts relevant information
#input: 
#   - single line from the kraken counts file
#returns:
#   - genome taxonomy ID
#   - total number of kmers for that genome
#   - dictionary: {classification id:number of kmers classified at this level}  
def parse_single_genome(curr_str):
    split_str = curr_str.strip().split('\t')
    #Error check for only one value in the line 
    if len(split_str) < 5:
        return [0,0,0]
    #Get the genome taxonomy ID tested
    genome_taxid = split_str[1]
    #Get which taxonomy ID this read length kmer mapped to
    mapped_id_kmers = {}
    total_kmers = 0
    kmer_distr = split_str[4]
    for kmers in kmer_distr.split():
        #Error check for no mapping/incorrect input file format
        if len(kmers.split(':')) == 1:
            #print("No kmers listed for this mapped_taxid: " + curr_str.strip())
            continue
        [curr_m_id, curr_kmers] = kmers.split(':')
        total_kmers += int(curr_kmers)
        if curr_m_id in mapped_id_kmers:
            mapped_id_kmers[curr_m_id] += int(curr_kmers)
        else:
            mapped_id_kmers[curr_m_id] = int(curr_kmers)
    #Error check for mappings
    if len(mapped_id_kmers) == 0:
        #print("No mappings for this taxid: " + curr_str.strip())
        return [0,0,0]
    #Return if correct
    return [genome_taxid, total_kmers, mapped_id_kmers]

#Main Method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', required=True,
        help='Kraken counts file for each genome mapped to the overall database.')
    parser.add_argument('-o', '--output', dest='output', required=True,
        help='Output file containing each classified taxonomy ID and the \
        kmer distributions of all genomes with this classification.')
    args=parser.parse_args()

    #Start Program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM START TIME: " + time + '\n')

    #Read in all input lines and create a dictionary mapping genome:mapped taxID:reads
    genome_dict = {}
    genome_dict_totalkmers = {}
    num_genomes = 0
    i_file = open(args.input, 'r')
    for line in i_file:
        [genome_taxid, total_kmers, mapped_taxids_kmers] = parse_single_genome(line)
        #No classification - ignore  
        if genome_taxid == 0:
            continue
        #Access genome in dictionary
        if genome_taxid not in genome_dict:
            genome_dict[genome_taxid] = {} 
            genome_dict_totalkmers[genome_taxid] = total_kmers
            num_genomes += 1
        else:
            genome_dict_totalkmers[genome_taxid] += total_kmers
        #Either save new kmer count or add to existing
        for m_taxid in mapped_taxids_kmers:
            if m_taxid not in genome_dict[genome_taxid]:
                genome_dict[genome_taxid][m_taxid] = mapped_taxids_kmers[m_taxid]
            else:
                genome_dict[genome_taxid][m_taxid] += mapped_taxids_kmers[m_taxid]
    i_file.close()
    sys.stdout.write('...' + str(num_genomes) + ' total genomes read from kraken output file\n')

    #Remap dictionary to create mapped_taxids:genome_taxid:kmers
    sys.stdout.write('...creating kmer counts file -- lists the number of kmers of each classification per genome\n')
    mapped_taxids_dict = {}
    for genome in genome_dict:
        for m_taxid in genome_dict[genome]:
            #Create new dictionary if needed
            if m_taxid not in mapped_taxids_dict:
                mapped_taxids_dict[m_taxid] = {}
            #Save number of kmers
            mapped_taxids_dict[m_taxid][genome] = genome_dict[genome][m_taxid]
    
    #Output distributions to file
    sys.stdout.write('...creating kmer distribution file -- lists genomes and kmer counts contributing to each genome\n')
    o_file = open(args.output, 'w')
    o_file.write('mapped_taxid\t' + 'genome_taxids:kmers_mapped:total_genome_kmers\n')
    for m_taxid in mapped_taxids_dict:
        o_file.write(m_taxid + '\t') 
        for genome_taxid in mapped_taxids_dict[m_taxid]:
            o_file.write(genome_taxid + ':' + str(mapped_taxids_dict[m_taxid][genome_taxid]))
            o_file.write(':' + str(genome_dict_totalkmers[genome_taxid]))
            o_file.write(' ')
        o_file.write('\n')
    o_file.close()

    #End of program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM END TIME: " + time + '\n')

if __name__ == "__main__":
    main()
