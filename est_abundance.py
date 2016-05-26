#!/usr/bin/python

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
import sys, argparse

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
    #Extract relevant information
    all_reads =  int(split_str[1])
    level_reads = int(split_str[2])
    level_type = split_str[3]
    taxid = split_str[4] 
    #Get name and spaces
    spaces = 0
    name = split_str[5]
    for char in name:
        if char == ' ':
            name = name[1:]
            spaces += 1 
        else:
            break 
    #Determine which level based on number of spaces
    level_num = spaces/2
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
    
    #Initialize variables
    root_node = -1
    prev_node = -1
    main_lvls = ['R','K','D','P','C','O','F','G','S']

    #Parse kraken report file and create tree 
    i_file = open(args.input, 'r')
    map2lvl_taxids = {}
    lvl_taxids = {} 
    last_taxid = -1
    for line in i_file:
        [name, taxid, level_num, level_id, all_reads, level_reads] = process_kraken_report(line)
        #Skip unclassified 
        if level_id == 'U':
            continue
        #Tree Root 
        if taxid == '1':
            root_node = Tree(name, taxid, level_num, 'R', all_reads, level_reads)
            prev_node = root_node
            continue 
        #Move to correct parent
        while level_num < (prev_node.level_num + 1):
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
            #Account for threshold at level
            if all_reads < int(args.thresh):
                last_taxid = -1
                continue
            #If level contains enough reads - save for abundance estimation
            lvl_taxids[taxid] = [name, all_reads, level_reads,0]
            last_taxid = taxid
            map2lvl_taxids[taxid] = [taxid, all_reads,0]
        elif main_lvls.index(level_id[0]) >= main_lvls.index(args.level):
            #For all nodes below the desired level 
            if last_taxid == -1:
                continue
            map2lvl_taxids[taxid] = [last_taxid, all_reads,0]
        #Add node to tree
        curr_node = Tree(name, taxid, level_num, level_id, all_reads, level_reads, None, prev_node)
        prev_node.add_child(curr_node)
        prev_node = curr_node 
    i_file.close()
    
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
    while len(curr_nodes) > 0:
        curr_node = curr_nodes.pop(0)
        #For each child node, if not at level, add to list of nodes to evaluate 
        for child_node in curr_node.children:
            if child_node.level_id != args.level:
                curr_nodes.append(child_node) 
        #Evaluate
        if curr_node.taxid not in kmer_distr_dict:
            continue
        #Get the dictionary listing all genomes mapping to this node
        curr_dict = kmer_distr_dict[curr_node.taxid]
        probability_dict_prelim = {}                  
        all_genome_reads = 0
        for genome in curr_dict:
            #Get the fraction of kmers of the genome expected to map to this node 
            fraction = float(curr_dict[genome][0])
            
            #Determine the number of reads classified by Kraken uniquely for the genome
            #and the fraction of the genome that is unique
            num_classified_reads = map2lvl_taxids[genome][1]
            if genome in kmer_distr_dict and genome in kmer_distr_dict[genome]:
                lvl_fraction = float(kmer_distr_dict[genome][genome][0])
            else:
                lvl_fraction = 1
            #Based on the classified reads and the fraction of unique reads, estimate
            #the true number of reads belonging to this genome in the sample 
            est_genome_reads = num_classified_reads/lvl_fraction
            all_genome_reads += est_genome_reads
            
            #Save values
            probability_dict_prelim[genome] = [fraction, est_genome_reads]
       
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

    #Print modified kraken report 
    #new_report, extension = os.path.splitext(args.input)
    #report_file = open(new_report + '_bracken' + extension, 'w')
    #report_file.close() 

if __name__ == "__main__":
    main()
