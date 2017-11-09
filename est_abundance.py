#!/usr/bin/python
#####################################################################
#est_abundance.py calculates abundances for a given microbial dataset classified by Kraken
#Copyright (C) 2016-2017 Jennifer Lu, jlu26@jhmi.edu

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
#02/29/2016
#
#This program takes the kraken report output and makes
#species or genus level abundance estimates using a Bayesian
#method relying on the Kraken assigned reads and
#the expected kmer distribution for each genome.
#
#Input files:
#   - Kraken report file generated from the kraken raw output file
#   - Kmer distribution file detailing the taxonomy IDs that each 
#     species/strain maps to when classified by kraken along with the
#     number of kmers out of total kmers mapping to the taxonomy ID
#
#Additional Input Parameters to Specify [OPTIONAL]:
#   - Classification Level to estimate abundances for
#     [Possible options = K, P, C, O, F, G, S, Default = S]  
#   - Read Threshold required to give confidence that higher 
#     classification level reads belong to a given species/strain
#     [Default = 10]
#
#Output file Format (tab-delimited, for species abundance estimation)
#   - Name of species
#   - Taxonomy ID for that species
#   - Taxonomy Level specified for (in default case 'S')
#   - Number of reads Kraken classified at that species
#   - Number of reads added by this abundance estimation
#   - Number of reads estimated for that species
#   - Fraction of total reads in the sample estimated for this species
#
#Methods:
#   - main
#   - process_kmer_distribution  
#   - process_kraken_report 
#
#####################################################################
import os, sys, argparse
import operator
from time import gmtime
from time import strftime

#Tree class
#usage: tree node used in constructing a taxonomy tree
#   including only the taxonomy levels and genomes identified in the Kraken report 
class Tree(object):
    'Tree node.'
    def __init__(self, name, taxid, level_num, level_id, all_reads, lvl_reads, children=None, parent=None):
        self.name = name
        self.taxid = taxid
        self.level_num = level_num
        self.level_id = level_id
        self.all_reads = all_reads
        self.lvl_reads = lvl_reads
        self.children = []
        self.parent =  parent
        if children is not None:
            for child in children:
                self.add_child(child)

    def add_child(self, node):
        assert isinstance(node,Tree)
        self.children.append(node) 

#process_kmer_distribution
#usage: parses a single line in the kmer distribution file and extracts
#relevant information for the genomes in this sample
#input:
#   - kmer distribution file generaed by generate_kmer_distribution.py.
#returns:
#   - classification taxonomy ID for this line
#   - dictionary of genomes/fractions of the genomes mapping to this classification
def process_kmer_distribution(curr_str, lvl_taxids, map2lvl_taxids):
    split_str = curr_str.strip().split('\t')
    #Parse each genome taxid mapping to m_taxid -- create dictionary instance
    temp_dict = {}
    mapped_taxid = split_str[0]
    for genome_str in split_str[1].split(' '):
        [g_taxid,mkmers,tkmers] = genome_str.split(':')
        mkmers = float(mkmers)
        tkmers = float(tkmers)
        fraction = mkmers/tkmers
        #Only include mappings for genomes within this sample
        if g_taxid in lvl_taxids or g_taxid in map2lvl_taxids:
            if g_taxid not in temp_dict:
                temp_dict[g_taxid] = [fraction]
            else:
                temp_dict[g_taxid].append(fraction)
    #Error check for relevant classifications 
    if len(temp_dict) == 0:
        return [-1,{}]
    #Return dictionary 
    return [mapped_taxid, temp_dict]

#process_kraken_report
#usage: parses a single line in the kraken report and extracts relevant information
#input: kraken report file with the following tab delimited lines 
#   - percent of total reads   
#   - number of reads (including at lower levels)
#   - number of reads (only at this level)
#   - taxonomy classification of level 
#       (U, - (root), - (cellular org), D, P, C, O, F, G, S) 
#   - taxonomy ID (0 = unclassified, 1 = root, 2 = Bacteria...etc)
#   - spaces + name 
#returns:
#   - classification/genome name
#   - taxonomy ID for this classification
#   - level for this classification (number)
#   - level name (U, -, D, P, C, O, F, G, S)
#   - all reads classified at this level and below in the tree
#   - reads classified only at this level
def process_kraken_report(curr_str):
    split_str = curr_str.strip().split('\t')
    try:
        int(split_str[1])
    except ValueError:
        return []
    #Extract relevant information
    all_reads =  int(split_str[1])
    level_reads = int(split_str[2])
    level_type = split_str[3]
    taxid = split_str[4] 
    #Get name and spaces
    spaces = 0
    name = split_str[-1]
    for char in name:
        if char == ' ':
            name = name[1:]
            spaces += 1 
        else:
            break 
    #Determine which level based on number of spaces
    level_num = int(spaces/2)
    return [name, taxid, level_num, level_type, all_reads, level_reads]
    
#Main method 
def main():
    #Parse arguments
    parser = argparse.ArgumentParser() 
    parser.add_argument('-i' ,'--input', dest='input', required=True,
        help='Input kraken report file.')
    parser.add_argument('-k', '--kmer_distr', dest='kmer_distr', required=True,
        help='Kmer distribution file.')
    parser.add_argument('-o', '--output', dest='output', required=True,
        help='Output modified kraken report file with abundance estimates')
    parser.add_argument('-l', '--level', dest='level', required=False,
        default='S',
        choices=['K','P','C','O','F','G','S'],
        help='Level to push all reads to.')
    parser.add_argument('-t', '--thresh','--threshold',dest='thresh', 
        required=False,default=10,
        help='Threshold for the minimum number of reads kraken must assign\
        to a classification for that classification to be considered in the\
        final abundance estimation.') 
    args=parser.parse_args()
    
    #Start program 
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM START TIME: " + time + '\n')

    #Abundance level
    lvl_dict = {}
    lvl_dict['K'] = 'kingdoms'
    lvl_dict['P'] = 'phylums'
    lvl_dict['O'] = 'orders'
    lvl_dict['C'] = 'classes'
    lvl_dict['F'] = 'families'
    lvl_dict['G'] = 'genuses' 
    lvl_dict['S'] = 'species'
    abundance_lvl = lvl_dict[args.level]

    #Initialize variables
    root_node = -1
    prev_node = -1
    main_lvls = ['R','K','D','P','C','O','F','G','S']
    total_reads = 0
    kept_reads = 0
    ignored_reads = 0
    n_lvl_total = 0
    n_lvl_est = 0
    n_lvl_del = 0
    lvl_nodes = []
    leaf_nodes = []

    #Parse kraken report file and create tree 
    i_file = open(args.input, 'r')
    map2lvl_taxids = {}
    lvl_taxids = {} 
    last_taxid = -1
    for line in i_file:
        report_vals = process_kraken_report(line)
        if len(report_vals) < 5:
            continue
        [name, taxid, level_num, level_id, all_reads, level_reads] = report_vals
        total_reads += level_reads
        #Skip unclassified 
        if level_id == 'U':
            unclassified_line = line
            u_reads = level_reads
            continue
        #Tree Root 
        if taxid == '1':
            root_node = Tree(name, taxid, level_num, 'R', all_reads, level_reads)
            prev_node = root_node
            continue 
        #Save leaf nodes
        if level_num != (prev_node.level_num + 1):
            leaf_nodes.append(prev_node)
        #Move to correct parent
        while level_num != (prev_node.level_num + 1):
            prev_node = prev_node.parent
        #Determine correct level ID
        if level_id == '-':
            if prev_node.level_id in main_lvls:
                level_id = prev_node.level_id + '1'
            else:
                num = int(prev_node.level_id[-1]) + 1
                level_id = prev_node.level_id[:-1] + str(num)
        #Desired level for abundance estimation or below
        if level_id == args.level:
            n_lvl_total += 1
            #Account for threshold at level
            if all_reads < int(args.thresh):
                n_lvl_del += 1
                ignored_reads += all_reads
                last_taxid = -1
            else: 
                #If level contains enough reads - save for abundance estimation
                n_lvl_est += 1
                kept_reads += all_reads
                lvl_taxids[taxid] = [name, all_reads, level_reads, 0]
                last_taxid = taxid
                map2lvl_taxids[taxid] = [taxid, all_reads, 0]
        elif main_lvls.index(level_id[0]) >= main_lvls.index(args.level):
            #For all nodes below the desired level 
            if last_taxid != -1:
                map2lvl_taxids[taxid] = [last_taxid, all_reads,0]
        #Add node to tree
        curr_node = Tree(name, taxid, level_num, level_id, all_reads, level_reads, None, prev_node)
        prev_node.add_child(curr_node)
        prev_node = curr_node 
    i_file.close()
    #Add last node
    leaf_nodes.append(prev_node)
    
    #Read in kmer distribution file
    k_file = open(args.kmer_distr,'r')
    kmer_distr_dict = {}
    for line in k_file.readlines()[1:]:
        [mapped_taxid, mapped_taxid_dict] = process_kmer_distribution(line,lvl_taxids,map2lvl_taxids)
        if len(mapped_taxid_dict) == 0:
            continue
        kmer_distr_dict[mapped_taxid] = mapped_taxid_dict
    k_file.close() 

    #For each current parent node, distribute level reads to genomes
    curr_nodes = [root_node]
    nondistributed_reads = 0
    distributed_reads = 0
    lvl_reads = 0
    while len(curr_nodes) > 0:
        curr_node = curr_nodes.pop(0)
        #For each child node, if not at level, add to list of nodes to evaluate 
        for child_node in curr_node.children:
            if child_node.level_id != args.level:
                curr_nodes.append(child_node) 
        #If no level taxids (or below) produce this classification 
        if curr_node.lvl_reads == 0:
            continue 
        if curr_node.taxid not in kmer_distr_dict:
            nondistributed_reads += curr_node.lvl_reads
            continue
        #Get the dictionary listing all genomes mapping to this node
        distributed_reads += curr_node.lvl_reads
        curr_dict = kmer_distr_dict[curr_node.taxid]
        probability_dict_prelim = {}                  
        all_genome_reads = 0
        for genome in curr_dict:
            #Get the fraction of kmers of the genome expected to map to this node 
            fraction = float(curr_dict[genome][0])
            
            #Determine the number of reads classified by Kraken uniquely for the genome
            #and the fraction of the genome that is unique
            num_classified_reads = float(map2lvl_taxids[genome][1])
            if genome in kmer_distr_dict and genome in kmer_distr_dict[genome]:
                lvl_fraction = float(kmer_distr_dict[genome][genome][0])
            else:
                lvl_fraction = 1.
            #Based on the classified reads and the fraction of unique reads, estimate
            #the true number of reads belonging to this genome in the sample 
            est_genome_reads = num_classified_reads/lvl_fraction
            all_genome_reads += est_genome_reads
            
            #Save values
            probability_dict_prelim[genome] = [fraction, est_genome_reads]
      
        if all_genome_reads == 0:
            continue
              
        #Get final probabilities
        #P_R_A = probability that a read is classified at the node given that it belongs to genome A
        #P_A = probability that a randomly selected read belongs to genome A
        #P_A_R = probability that a read belongs to genome A given that its classified at the node 
        total_probability = 0.0
        probability_dict_final = {}
        for genome in probability_dict_prelim:
            [P_R_A, est_g_reads] = probability_dict_prelim[genome]
            P_A = float(est_g_reads)/float(all_genome_reads)
            P_A_R = float(P_R_A)*float(P_A)
            probability_dict_final[genome] = P_A_R
            total_probability += P_A_R

        #Find the normalize probabilty and Distribute reads accordingly
        for genome in probability_dict_final:
            add_fraction = probability_dict_final[genome]/total_probability
            add_reads = add_fraction*float(curr_node.lvl_reads)
            map2lvl_taxids[genome][2] += add_reads
    
    #For all genomes, map reads up to level
    for genome in map2lvl_taxids:
        [lvl_taxid,all_reads,add_reads] = map2lvl_taxids[genome]
        lvl_taxids[lvl_taxid][3] += add_reads 

    #Sum all of the reads for the desired level -- use for fraction of reads
    sum_all_reads = 0
    for taxid in lvl_taxids:
        [name, all_reads, lvl_reads, added_reads] = lvl_taxids[taxid]
        new_all_reads = float(all_reads) + float(added_reads)
        sum_all_reads += new_all_reads

    #Print for each classification level: 
    #   - name, taxonomy ID, taxonomy level
    #   - kraken assigned reads, added reads, estimated reads, and fraction total reads 
    o_file = open(args.output, 'w')
    o_file.write('name\t' + 'taxonomy_id\t' + 'taxonomy_lvl\t' + 'kraken_assigned_reads\t' + 'added_reads\t' + 'new_est_reads\t' + 'fraction_total_reads\n')
    for taxid in lvl_taxids:
        [name, all_reads, lvl_reads, added_reads] = lvl_taxids[taxid]
        #Count up all added reads + all_reads already at the level
        new_all_reads = float(all_reads) + float(added_reads)
        #Output
        o_file.write(name + '\t')
        o_file.write(taxid + '\t')
        o_file.write(args.level + '\t')
        o_file.write(str(int(all_reads)) + '\t')
        o_file.write(str(int(new_all_reads)-int(all_reads))+'\t')
        o_file.write(str(int(new_all_reads)) + '\t')
        o_file.write("%0.5f\n" % (float(new_all_reads)/float(sum_all_reads)))
    o_file.close()
    
    #Print to screen
    print("BRACKEN SUMMARY (Kraken report: %s)" % args.input)
    print("    >>> Threshold: %i " % int(args.thresh))
    print("    >>> Number of %s in sample: %i " % (abundance_lvl, n_lvl_total))
    print("\t  >> Number of %s with reads > threshold: %i " % (abundance_lvl, n_lvl_est))
    print("\t  >> Number of %s with reads < threshold: %i " % (abundance_lvl, n_lvl_del))
    print("    >>> Total reads in sample: %i" % total_reads)
    print("\t  >> Total reads kept at %s level (reads > threshold): %i" %(abundance_lvl, kept_reads))
    print("\t  >> Total reads discarded (%s reads < threshold): %i" % (abundance_lvl, ignored_reads))
    print("\t  >> Reads distributed: %i" % distributed_reads)
    print("\t  >> Reads not distributed (eg. no %s above threshold): %i" % (abundance_lvl, nondistributed_reads))
    print("\t  >> Unclassified reads: %i" % u_reads)
    print("BRACKEN OUTPUT PRODUCED: %s" % args.output)
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM END TIME: " + time + '\n')
    
    ###########################################################################
    #Kraken-Style Report Section added 05/26/2016
    #Jennifer Lu, jlu26 
    #For each child node, add reads to all parents 
    new_reads = {}
    test_reads = 0
    for curr_leaf in leaf_nodes:
        #Move to estimation level
        curr_node = curr_leaf
        if args.level in curr_node.level_id:
            while args.level != curr_node.level_id:
                curr_node = curr_node.parent
        #Determine number of reads to add
        add_reads = curr_node.all_reads
        if curr_node.taxid in lvl_taxids:
            [name, all_reads, lvl_reads, added_reads] = lvl_taxids[curr_node.taxid]
            add_reads += added_reads
        #If this level tree already traversed, do not traverse
        if curr_node.taxid in new_reads:
            continue
        #Save reads for this node
        new_reads[curr_node.taxid] = add_reads
        test_reads += add_reads 
        #Traverse tree 
        while curr_node.parent is not None:
            #Move to parent
            curr_node = curr_node.parent
            #Add to dictionary if not previously found
            if curr_node.taxid not in new_reads:
                if curr_node.taxid not in kmer_distr_dict:
                    add_reads += curr_node.lvl_reads
                new_reads[curr_node.taxid] = 0
            #Add reads
            new_reads[curr_node.taxid] += add_reads 
    #Print modified kraken report 
    new_report, extension = os.path.splitext(args.input)
    r_file = open(new_report + '_bracken' + extension, 'w')
    r_file.write(unclassified_line)
    #For each current parent node, print to file 
    curr_nodes = [root_node]
    while len(curr_nodes) > 0:
        curr_node = curr_nodes.pop(0)
        #For each child node, add to list of nodes to evaluate 
        children = 0
        for child_node in sorted(curr_node.children, key=operator.attrgetter('all_reads')):
            #Add if at level or above 
            if child_node.level_id[0] != args.level or child_node.level_id == args.level:
                curr_nodes.insert(0,child_node) 
                children += 1
        #Print information for this level 
        #For level where estimate is made
        if curr_node.taxid in lvl_taxids:
            [name, all_reads, lvl_reads, added_reads] = lvl_taxids[curr_node.taxid]
            new_all_reads = float(all_reads) + float(added_reads)
       
        #Print information for this level
        new_all_reads = new_reads[curr_node.taxid]
        r_file.write("%0.2f\t" % (float(new_all_reads)/float(total_reads)*100))
        r_file.write("%i\t" % (new_all_reads))
        if children == 0:
            r_file.write("%i\t" % (new_all_reads))
        else:
            r_file.write("0\t")
        r_file.write(curr_node.level_id + "\t")
        r_file.write(curr_node.taxid + "\t")
        r_file.write(" "*curr_node.level_num*2 + curr_node.name + "\n")
    r_file.close() 
    ###########################################################################

if __name__ == "__main__":
    main()
